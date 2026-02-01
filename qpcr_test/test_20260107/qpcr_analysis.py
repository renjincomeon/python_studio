import json
import pprint
from statistics import mean, stdev
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import os

def calculate_delta_ct_v2(df_data, dict_info):
    """使用更简洁的方式计算ΔCt"""
    delta_ct_dict = {}
    
    for sample, genes in dict_info.items():
        sample_data = {}
        
        # 获取Actin平均值
        actin_positions = genes.get('Actin', [])
        if not actin_positions:
            print(f"跳过样本 {sample}：缺少Actin数据")
            continue
            
        try:
            actin_mean = df_data.loc[actin_positions, 'Cp'].mean()
        except KeyError as e:
            print(f"样本 {sample} 的Actin位置错误：{e}")
            continue
        
        # 计算每个基因的ΔCt
        for gene, positions in genes.items():
            if gene == 'Actin':
                continue
                
            if positions:
                try:
                    ct_values = df_data.loc[positions, 'Cp']
                    delta_ct_values = (ct_values - actin_mean).round(2).tolist()
                    sample_data[gene] = delta_ct_values
                except KeyError as e:
                    print(f"样本 {sample} 基因 {gene} 位置错误：{e}")
        
        if sample_data:
            delta_ct_dict[sample] = sample_data
    
    return delta_ct_dict

def calculate_delta_delta_ct_v2(delta_ct_dict, control_sample):
    """计算ΔΔCt"""
    if control_sample not in delta_ct_dict:
        print(f"错误：控制样本 {control_sample} 不存在于数据中")
        return None
    
    # 计算对照样本的ΔCt均值
    control_means = {}
    for gene, values in delta_ct_dict[control_sample].items():
        if values:  # 确保有数据
            control_means[gene] = round(mean(values), 2)
    
    # 计算所有样本的ΔΔCt
    delta_delta_ct_dict = {}
    for sample, genes in delta_ct_dict.items():
        sample_results = {}
        for gene, delta_ct_values in genes.items():
            if gene in control_means and delta_ct_values:
                delta_delta_ct_values = [
                    round(delta_ct - control_means[gene], 2) 
                    for delta_ct in delta_ct_values
                ]
                sample_results[gene] = delta_delta_ct_values
        
        if sample_results:
            delta_delta_ct_dict[sample] = sample_results
    
    return delta_delta_ct_dict

def get_fc(delta_delta_ct_dict):
    """计算Fold Change"""
    fc_dict = {}
    for sample, genes in delta_delta_ct_dict.items():
        for gene, values in genes.items():
            fold_changes = [round(2 ** (-val), 2) for val in values]
            fc_dict.setdefault(sample, {})[gene] = fold_changes
    return fc_dict

def create_summary_dataframe(delta_ct_dict, delta_delta_ct_dict, fc_dict):
    """创建汇总的DataFrame"""
    records = []
    
    for sample, genes in delta_delta_ct_dict.items():
        for gene, values in genes.items():
            # 获取ΔCt数据
            delta_ct_vals = delta_ct_dict.get(sample, {}).get(gene, [])
            delta_ct_mean = round(mean(delta_ct_vals), 2) if delta_ct_vals else None
            delta_ct_sem = round(stdev(delta_ct_vals)/np.sqrt(len(delta_ct_vals)), 2) if len(delta_ct_vals) > 1 else None
            
            # 获取ΔΔCt数据
            delta_delta_ct_mean = round(mean(values), 2) if values else None
            delta_delta_ct_sem = round(stdev(values)/np.sqrt(len(values)), 2) if len(values) > 1 else None
            
            # 获取Fold Change数据
            fc_vals = fc_dict.get(sample, {}).get(gene, [])
            fc_mean = round(mean(fc_vals), 2) if fc_vals else None
            fc_sem = round(stdev(fc_vals)/np.sqrt(len(fc_vals)), 2) if len(fc_vals) > 1 else None
            
            records.append({
                'Sample': sample,
                'Gene': gene,
                'DeltaCt_Mean': delta_ct_mean,
                'DeltaCt_SEM': delta_ct_sem,
                'DeltaDeltaCt_Mean': delta_delta_ct_mean,
                'DeltaDeltaCt_SEM': delta_delta_ct_sem,
                'FoldChange_Mean': fc_mean,
                'FoldChange_SEM': fc_sem,
                'Replicates': len(values)
            })
    
    return pd.DataFrame(records)

def plot_fold_change_bar(df_summary, output_dir="results/plots"):
    """绘制Fold Change柱状图"""
    os.makedirs(output_dir, exist_ok=True)
    
    fig, ax = plt.subplots(figsize=(12, 6))
    
    # 准备数据
    plot_df = df_summary.pivot(index='Sample', columns='Gene', values='FoldChange_Mean')
    error_df = df_summary.pivot(index='Sample', columns='Gene', values='FoldChange_SEM')
    
    # 绘制柱状图
    x = np.arange(len(plot_df.index))
    width = 0.8 / len(plot_df.columns)
    
    for i, gene in enumerate(plot_df.columns):
        offset = (i - len(plot_df.columns)/2 + 0.5) * width
        bars = ax.bar(x + offset, plot_df[gene].values, width, 
                     label=gene, alpha=0.8)
        
        # 添加误差棒
        if error_df is not None:
            ax.errorbar(x + offset, plot_df[gene].values, 
                       yerr=error_df[gene].values, 
                       fmt='none', color='black', capsize=3)
    
    ax.set_xlabel('Samples', fontsize=12)
    ax.set_ylabel('Fold Change (2^(-ΔΔCt))', fontsize=12)
    ax.set_title('qPCR Results - Gene Expression Changes', fontsize=14, fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels(plot_df.index, rotation=45, ha='right')
    ax.legend(title='Genes', bbox_to_anchor=(1.05, 1), loc='upper left')
    ax.axhline(y=1, color='gray', linestyle='--', alpha=0.5)
    ax.set_yscale('log')
    
    plt.tight_layout()
    plt.savefig(f'{output_dir}/fold_change_barplot.png', dpi=300, bbox_inches='tight')
    plt.savefig(f'{output_dir}/fold_change_barplot.pdf', bbox_inches='tight')
    plt.show()

def plot_gene_expression_heatmap(df_summary, output_dir="results/plots"):
    """绘制基因表达热图"""
    os.makedirs(output_dir, exist_ok=True)
    
    # 准备数据
    pivot_df = df_summary.pivot(index='Sample', columns='Gene', values='FoldChange_Mean')
    
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # 使用seaborn绘制热图
    sns.heatmap(pivot_df, 
                annot=True, 
                fmt='.2f',
                cmap='RdBu_r',
                center=1,
                linewidths=0.5,
                cbar_kws={'label': 'Fold Change'},
                ax=ax)
    
    ax.set_title('Gene Expression Heatmap (Fold Change)', fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(f'{output_dir}/gene_expression_heatmap.png', dpi=300, bbox_inches='tight')
    plt.savefig(f'{output_dir}/gene_expression_heatmap.pdf', bbox_inches='tight')
    plt.show()

def plot_individual_replicates(delta_delta_ct_dict, fc_dict, output_dir="results/plots"):
    """绘制每个样本的重复数据点"""
    os.makedirs(output_dir, exist_ok=True)
    
    # 准备数据
    replicate_data = []
    for sample, genes in delta_delta_ct_dict.items():
        for gene, values in genes.items():
            for i, (ddct, fc) in enumerate(zip(values, fc_dict.get(sample, {}).get(gene, []))):
                replicate_data.append({
                    'Sample': sample,
                    'Gene': gene,
                    'Replicate': i+1,
                    'DeltaDeltaCt': ddct,
                    'FoldChange': fc
                })
    
    df_replicates = pd.DataFrame(replicate_data)
    
    # 绘制散点图
    fig, axes = plt.subplots(1, 2, figsize=(15, 6))
    
    # ΔΔCt散点图
    ax1 = axes[0]
    for gene in df_replicates['Gene'].unique():
        gene_data = df_replicates[df_replicates['Gene'] == gene]
        ax1.scatter(gene_data['Sample'], gene_data['DeltaDeltaCt'], 
                   alpha=0.6, s=100, label=gene)
    
    ax1.set_xlabel('Samples', fontsize=12)
    ax1.set_ylabel('ΔΔCt', fontsize=12)
    ax1.set_title('ΔΔCt Replicate Data Points', fontsize=14, fontweight='bold')
    ax1.legend(title='Genes', bbox_to_anchor=(1.05, 1), loc='upper left')
    ax1.tick_params(axis='x', rotation=45)
    ax1.axhline(y=0, color='gray', linestyle='--', alpha=0.5)
    
    # Fold Change散点图
    ax2 = axes[1]
    for gene in df_replicates['Gene'].unique():
        gene_data = df_replicates[df_replicates['Gene'] == gene]
        ax2.scatter(gene_data['Sample'], gene_data['FoldChange'], 
                   alpha=0.6, s=100, label=gene)

    ax2.set_xlabel('Samples', fontsize=12)
    ax2.set_ylabel('Fold Change', fontsize=12)
    ax2.set_title('Fold Change Replicate Data Points', fontsize=14, fontweight='bold')
    ax2.legend(title='Genes', bbox_to_anchor=(1.05, 1), loc='upper left')
    ax2.tick_params(axis='x', rotation=45)
    ax2.set_yscale('log')
    ax2.axhline(y=1, color='gray', linestyle='--', alpha=0.5)
    
    plt.tight_layout()
    plt.savefig(f'{output_dir}/replicates_scatterplot.png', dpi=300, bbox_inches='tight')
    plt.savefig(f'{output_dir}/replicates_scatterplot.pdf', bbox_inches='tight')
    plt.show()

def plot_statistical_comparison(df_summary, control_sample='sample4', output_dir="results/plots"):
    """绘制统计比较图"""
    os.makedirs(output_dir, exist_ok=True)
    
    genes = df_summary['Gene'].unique()
    n_genes = len(genes)
    
    fig, axes = plt.subplots(n_genes, 1, figsize=(10, 4*n_genes))
    if n_genes == 1:
        axes = [axes]
    
    for idx, gene in enumerate(genes):
        ax = axes[idx]
        gene_data = df_summary[df_summary['Gene'] == gene]
        
        # 获取对照组数据
        control_data = gene_data[gene_data['Sample'] == control_sample]
        control_fc = control_data['FoldChange_Mean'].values[0]
        
        # 绘制条形图
        samples = gene_data['Sample'].tolist()
        means = gene_data['FoldChange_Mean'].tolist()
        sems = gene_data['FoldChange_SEM'].tolist()
        
        bars = ax.bar(range(len(samples)), means, yerr=sems, 
                     capsize=5, alpha=0.7, color='skyblue')
        
        # 标记对照组
        if control_sample in samples:
            control_idx = samples.index(control_sample)
            bars[control_idx].set_color('gray')
            bars[control_idx].set_alpha(0.5)
        
        # 添加显著性标记（简化版本，实际可能需要t-test）
        for i, (sample, mean_val) in enumerate(zip(samples, means)):
            if sample != control_sample:
                # 这里简化处理，实际应用中应进行统计检验
                diff = mean_val - control_fc
                if abs(diff) > control_fc * 0.5:  # 示例条件
                    ax.text(i, mean_val + sems[i] + 0.1, '*', 
                           ha='center', va='bottom', fontsize=14)
        
        ax.set_xticks(range(len(samples)))
        ax.set_xticklabels(samples, rotation=45, ha='right')
        ax.set_ylabel('Fold Change', fontsize=12)
        ax.set_title(f'{gene} Expression Level Comparison', fontsize=13, fontweight='bold')
        ax.axhline(y=1, color='gray', linestyle='--', alpha=0.5)
        ax.set_yscale('log')
    
    plt.tight_layout()
    plt.savefig(f'{output_dir}/statistical_comparison.png', dpi=300, bbox_inches='tight')
    plt.savefig(f'{output_dir}/statistical_comparison.pdf', bbox_inches='tight')
    plt.show()

def plot_quality_control(df_raw, dict_info, output_dir="results/plots"):
    """绘制质控图"""
    os.makedirs(output_dir, exist_ok=True)
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    axes = axes.flatten()
    
    # 1. Actin Ct值分布
    actin_cts = []
    labels = []
    for sample, genes in dict_info.items():
        actin_positions = genes.get('Actin', [])
        if actin_positions:
            ct_values = df_raw.loc[actin_positions, 'Cp'].tolist()
            actin_cts.extend(ct_values)
            labels.extend([sample] * len(ct_values))
    
    if actin_cts:
        ax1 = axes[0]
        sample_order = sorted(set(labels))
        positions = range(len(sample_order))
        
        # 创建箱线图
        box_data = []
        for sample in sample_order:
            sample_cts = [ct for ct, lab in zip(actin_cts, labels) if lab == sample]
            box_data.append(sample_cts)
        
        bp = ax1.boxplot(box_data, positions=positions, patch_artist=True)
        
        # 设置颜色
        colors = plt.cm.Set3(np.linspace(0, 1, len(sample_order)))
        for patch, color in zip(bp['boxes'], colors):
            patch.set_facecolor(color)
        
        ax1.set_xticks(positions)
        ax1.set_xticklabels(sample_order, rotation=45, ha='right')
        ax1.set_ylabel('Ct Value', fontsize=12)
        ax1.set_title('Actin Gene Ct value distribution (Internal Control Gene)', fontsize=13, fontweight='bold')
        ax1.grid(True, alpha=0.3)
    
    # 2. 技术重复一致性
    ax2 = axes[1]
    all_replicates = []
    for sample, genes in dict_info.items():
        for gene, positions in genes.items():
            if gene != 'Actin' and positions:
                ct_values = df_raw.loc[positions, 'Cp'].tolist()
                if len(ct_values) > 1:
                    cv = stdev(ct_values) / mean(ct_values) * 100 if mean(ct_values) != 0 else 0
                    all_replicates.append({
                        'Sample': sample,
                        'Gene': gene,
                        'CV': cv
                    })
    
    if all_replicates:
        rep_df = pd.DataFrame(all_replicates)
        sns.boxplot(data=rep_df, x='Gene', y='CV', ax=ax2)
        ax2.set_xlabel('Gene', fontsize=12)
        ax2.set_ylabel('Coefficient of Variation (%)', fontsize=12)
        ax2.set_title('Technical Replicate Consistency (CV%)', fontsize=13, fontweight='bold')
        ax2.tick_params(axis='x', rotation=45)
        ax2.axhline(y=10, color='r', linestyle='--', alpha=0.5, label='10% threshold')
        ax2.legend()
    
    # 3. 总体Ct值分布
    ax3 = axes[2]
    all_ct_values = []
    for sample, genes in dict_info.items():
        for gene, positions in genes.items():
            if positions:
                ct_values = df_raw.loc[positions, 'Cp'].tolist()
                all_ct_values.extend(ct_values)
    
    if all_ct_values:
        ax3.hist(all_ct_values, bins=20, alpha=0.7, edgecolor='black')
        ax3.set_xlabel('Ct Value', fontsize=12)
        ax3.set_ylabel('Frequency', fontsize=12)
        ax3.set_title('Overall Ct Value Distribution', fontsize=13, fontweight='bold')
        ax3.axvline(x=mean(all_ct_values), color='r', linestyle='--', 
                   label=f'Mean: {mean(all_ct_values):.2f}')
        ax3.legend()
    
    # 4. 样本间比较
    ax4 = axes[3]
    sample_means = []
    sample_names = []
    for sample, genes in dict_info.items():
        ct_values = []
        for gene, positions in genes.items():
            if positions:
                ct_values.extend(df_raw.loc[positions, 'Cp'].tolist())
        if ct_values:
            sample_means.append(mean(ct_values))
            sample_names.append(sample)
    
    if sample_means:
        ax4.bar(range(len(sample_names)), sample_means, alpha=0.7)
        ax4.set_xticks(range(len(sample_names)))
        ax4.set_xticklabels(sample_names, rotation=45, ha='right')
        ax4.set_ylabel('Average Ct Value', fontsize=12)
        ax4.set_title('Average Ct Value Comparison Between Samples', fontsize=13, fontweight='bold')
    
    plt.tight_layout()
    plt.savefig(f'{output_dir}/quality_control_plots.png', dpi=300, bbox_inches='tight')
    plt.savefig(f'{output_dir}/quality_control_plots.pdf', bbox_inches='tight')
    plt.show()

def save_results(delta_ct_dict, delta_delta_ct_dict, fc_dict, df_summary, output_dir="results"):
    """保存所有结果"""
    os.makedirs(output_dir, exist_ok=True)
    
    # 保存数据到JSON
    results = {
        'delta_ct': delta_ct_dict,
        'delta_delta_ct': delta_delta_ct_dict,
        'fold_change': fc_dict
    }
    
    with open(f'{output_dir}/qPCR_results.json', 'w', encoding='utf-8') as f:
        json.dump(results, f, indent=2, ensure_ascii=False)
    
    # 保存汇总表格到CSV和Excel
    df_summary.to_csv(f'{output_dir}/qPCR_summary.csv', index=False, encoding='utf-8-sig')
    
    with pd.ExcelWriter(f'{output_dir}/qPCR_analysis_results.xlsx') as writer:
        df_summary.to_excel(writer, sheet_name='Summary', index=False)
        
        # 添加详细数据表
        detailed_data = []
        for sample, genes in delta_ct_dict.items():
            for gene, values in genes.items():
                for i, (dct, ddct, fc) in enumerate(zip(
                    values,
                    delta_delta_ct_dict.get(sample, {}).get(gene, []),
                    fc_dict.get(sample, {}).get(gene, [])
                )):
                    detailed_data.append({
                        'Sample': sample,
                        'Gene': gene,
                        'Replicate': i+1,
                        'DeltaCt': dct,
                        'DeltaDeltaCt': ddct,
                        'FoldChange': fc
                    })
        
        df_detailed = pd.DataFrame(detailed_data)
        df_detailed.to_excel(writer, sheet_name='Detailed_Data', index=False)
    
    print(f"\n结果已保存到 {output_dir}/ 目录")
    print(f"- qPCR_results.json: 原始计算结果")
    print(f"- qPCR_summary.csv: 汇总表格")
    print(f"- qPCR_analysis_results.xlsx: Excel格式结果")
    print(f"- plots/: 各种统计图表")

def export_fc_to_excel(fc_dict, filename):
    dict_new = {}
    for key, value in fc_dict.items():
        for gene, fc_values in value.items():
            for i, fc in enumerate(fc_values):
                key_str = '{}_{}'.format(key, i+1)
                if key_str in dict_new:
                    dict_new[key_str].update({gene: fc})
                else:
                    dict_new[key_str] = {gene: fc}
    fc_df = pd.DataFrame(dict_new)
    fc_df.to_excel(filename)
    return fc_df


def main(df_raw, json_data, control_sample='sample4'):
    """
    主分析函数
    """
    print("=" * 60)
    print("qPCR数据分析开始")
    print("=" * 60)
    
    # 1. 计算ΔCt
    print("\n1. 计算ΔCt值...")
    delta_ct_dict = calculate_delta_ct_v2(df_raw, json_data)
    
    print(f"ΔCt计算结果：")
    for sample, genes in delta_ct_dict.items():
        print(f"  样本 {sample}:")
        for gene, values in genes.items():
            print(f"    {gene}: {values} (均值: {round(mean(values), 2)})")
    
    # 2. 计算ΔΔCt
    print(f"\n2. 计算ΔΔCt值（对照样本: {control_sample}）...")
    delta_delta_ct_dict = calculate_delta_delta_ct_v2(delta_ct_dict, control_sample)
    if delta_delta_ct_dict:
        print(f"ΔΔCt计算结果：")
        for sample, genes in delta_delta_ct_dict.items():
            print(f"  样本 {sample}:")
            for gene, values in genes.items():
                mean_val = round(mean(values), 2) if values else None
                print(f"    {gene}: {values} (均值: {mean_val})")
    
    # 3. 计算2^(-ΔΔCt)
    print(f"\n3. 计算2^(-ΔΔCt)值...")
    fc_dict = get_fc(delta_delta_ct_dict)
    export_fc_to_excel(fc_dict, filename="results/fold_change_values.xlsx")
    print(f"Fold Change计算结果：")
    for sample, genes in fc_dict.items():
        print(f"  样本 {sample}:")
        for gene, values in genes.items():
            mean_val = round(mean(values), 2) if values else None
            print(f"    {gene}: {values} (均值: {mean_val})")
    
    # 4. 创建汇总DataFrame
    print(f"\n4. 创建数据汇总...")
    df_summary = create_summary_dataframe(delta_ct_dict, delta_delta_ct_dict, fc_dict)
    print(f"数据汇总完成，共 {len(df_summary)} 条记录")
    
    # 5. 绘制统计图表
    print(f"\n5. 生成统计图表...")
    
    # 5.1 质控图
    print("   - 生成质控图...")
    plot_quality_control(df_raw, json_data)
    
    # 5.2 主要结果图
    print("   - 生成Fold Change柱状图...")
    plot_fold_change_bar(df_summary)
    
    print("   - 生成基因表达热图...")
    plot_gene_expression_heatmap(df_summary)
    
    print("   - 生成重复数据点图...")
    plot_individual_replicates(delta_delta_ct_dict, fc_dict)
    
    print("   - 生成统计比较图...")
    plot_statistical_comparison(df_summary, control_sample)
    
    # 6. 保存结果
    print(f"\n6. 保存分析结果...")
    save_results(delta_ct_dict, delta_delta_ct_dict, fc_dict, df_summary)
    
    print("\n" + "=" * 60)
    print("分析完成！")
    print("所有图表已保存至 results/plots/ 目录")
    print("所有数据已保存至 results/ 目录")
    print("=" * 60)
    
    return delta_ct_dict, delta_delta_ct_dict, fc_dict, df_summary


if __name__ == '__main__':
    # 创建输出目录
    wd_dir = os.path.dirname(sys.argv[0])
    os.chdir(wd_dir)
    os.makedirs('results/plots', exist_ok=True)
    
    # 加载数据
    try:
        with open('./config.json', 'r') as f:
            json_data = json.load(f)
        print("✓ 成功加载配置文件")
    except FileNotFoundError:
        print("错误：找不到config.json文件")
        exit(1)
    except json.JSONDecodeError:
        print("错误：配置文件格式不正确")
        exit(1)
    
    try:
        df_raw = pd.read_csv('./renjin_qpcr_20260107_2.txt', sep='\t', skiprows=1)
        df_raw.index = df_raw['Pos']
        print("✓ 成功加载qPCR数据")
        print(f"数据形状: {df_raw.shape}")
    except FileNotFoundError:
        print("错误：找不到数据文件")
        exit(1)
    except Exception as e:
        print(f"错误：加载数据时出错 - {e}")
        exit(1)
    
    # 运行分析
    try:
        results = main(df_raw, json_data, control_sample='EV+F-Ctrl')
        print("\n✓ 分析成功完成！")
    except Exception as e:
        print(f"\n✗ 分析过程中出错: {e}")
        import traceback
        traceback.print_exc()