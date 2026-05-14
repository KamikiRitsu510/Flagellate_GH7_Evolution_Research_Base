# treecolorize — 实战示例

> **场景**：你刚跑完 IQ-TREE，手里有一棵包含 534 条序列的 GH7 纤维素酶系统发育树。序列来自细菌、真菌、白蚁肠道鞭毛虫共生体等多个类群，tip label 格式是 `accession 物种名`（NCBI 标准格式）。你希望按生物学分类着色，生成一张可发表的 PDF 图，同时导出一个可编辑的 ggplot 对象用于后续调整。

---

## 0. 准备工作

```r
# 安装（首次）
install.packages("path/to/treecolorize", repos = NULL, type = "source")
install.packages("ape")
BiocManager::install(c("ggtree", "treeio"))   # 可选，ggtree 后端需要

library(ape)
library(treecolorize)
```

查看快速帮助：

```r
treecolorize_help()
```

输出示例：

```
── treecolorize 快速参考 ──────────────────────────────────────────────
主函数:
  treecolorize(tree, keyword_map, backend, output_file, ...)
  → 一站式：分类 → 配色 → 画图

分步函数:
  classify_tips(tip.labels, keyword_map)   → 分组向量
  auto_palette(groups)                     → 配色列表
  plot_tree_ape(tree, groups, palette)     → ape 出图
  plot_tree_ggtree(tree, groups, palette)  → 返回 ggplot 对象
  save_tree_ape(tree, groups, file)        → 保存 PDF/PNG

branch_method 选项: "majority" | "tip_only" | "uniform"
输入 ?treecolorize 查看完整文档。
──────────────────────────────────────────────────────────────────────
```

---

## 1. 读取树文件

```r
tree <- read.tree("GH7_final_534.treefile")

# 看看 tip label 格式
head(tree$tip.label, 6)
# [1] "BAB64565.3 Pseudotrichonympha_grassii"
# [2] "AAY83390.3 Pseudotrichonympha_sp."
# [3] "WP_165634028.1 Pseudobacteroides_cellulosivorans"
# [4] "XP_006965674.1 Reticulitermes_flavipes"
# [5] "CAD48174.1 Piromyces_sp."
# [6] "OXB12345.1 Nasutitermes_corniger"

cat("共", length(tree$tip.label), "条序列\n")
# 共 534 条序列
```

---

## 2. 定义分组规则

```r
# 关键词映射：在 tip label 中搜索关键词 → 分配组名
# 支持正则表达式，大小写不敏感
keyword_map <- list(
  # 鞭毛虫共生体（核心研究对象）
  "Pseudotrichonympha|Holomastigotoides|Spirotrichonympha|Reticulomyxa" = "Flagellate symbiont",
  
  # 白蚁宿主
  "Reticulitermes|Nasutitermes|Coptotermes|Macrotermes"                 = "Termite host",
  
  # 细菌来源（HGT 供体候选）
  "Pseudobacteroides|Ruminococcus|Fibrobacter|Clostridium|WP_|NZ_"     = "Bacteria",
  
  # 真菌（外类群）
  "Piromyces|Neocallimastix|Orpinomyces|Fungi|fungus"                  = "Fungus (outgroup)",
  
  # 其他真核生物
  "XP_|NP_"                                                             = "Other eukaryote"
)
```

> **提示**：`WP_` 和 `NZ_` 是 NCBI 细菌蛋白/基因组 accession 前缀，可以用作快速分类标志。`XP_`/`NP_` 是真核预测/精选蛋白前缀。

---

## 3. 分类并检查

```r
groups <- classify_tips(tree$tip.label, keyword_map)

# 检查分组结果
table(groups)
# groups
#          Bacteria Flagellate symbiont       Fungus (outgroup) 
#               287                  17                       8 
#     Other eukaryote       Termite host              Other 
#               198                  17                       7

# 检查"Other"组（未匹配的）——可能需要补充关键词
other_tips <- names(groups[groups == "Other"])
head(other_tips)
# [1] "CAB87654.2 uncultured_bacterium"
# → 补充关键词 "uncultured" = "Bacteria"
```

---

## 4. 配色

```r
palette <- auto_palette(
  groups,
  # 可选：手动指定某些组的颜色
  user_colors = c(
    "Flagellate symbiont" = "#E63946",   # 红色突出核心序列
    "Fungus (outgroup)"  = "#6A4C93"    # 紫色标外类群
  ),
  other_color = "#CCCCCC"
)

# 查看配色方案
print(palette)
# $`Flagellate symbiont`
# [1] "#E63946"
# $Bacteria
# [1] "#F4A261"
# ...
```

---

## 5. 用 ape 绘图并保存

```r
# 保存高质量 PDF（自动设置画布）
save_tree_ape(
  tree, groups,
  file    = "GH7_tree_colored.pdf",
  palette = palette,
  width   = 14,
  height  = 36,             # 序列多时适当加高
  type    = "phylogram",
  branch_method = "majority",   # 内部节点颜色 = 子节点多数投票
  tip_cex = 0.35,               # 序列多时缩小字号
  edge_width = 0.6,
  title   = "GH7 Cellulase Phylogeny (534 sequences)"
)

# 同时看预览
plot_tree_ape(tree, groups, palette,
              type = "fan",        # 扇形树，序列多时更紧凑
              tip_cex = 0.3,
              legend_pos = "bottomleft")
```

输出：`GH7_tree_colored.pdf`，各分支按生物学类群着色，图例自动生成。

---

## 6. 用 ggtree 后端生成可编辑对象

```r
library(ggtree)
library(ggplot2)

p <- plot_tree_ggtree(
  tree, groups, palette,
  layout          = "rectangular",
  tip_label_size  = 1.0,
  tip_point_size  = 0.6,
  label_offset    = 0.05,
  title           = "GH7 Cellulase Phylogeny"
)

# ggplot 对象可以继续叠加图层
p <- p +
  geom_treescale(x = 0, y = -5, fontsize = 3) +   # 添加比例尺
  theme(legend.text = element_text(size = 8))

ggsave("GH7_tree_ggtree.pdf", p, width = 14, height = 36)
```

---

## 7. 最终结果

```
输出文件：
  GH7_tree_colored.pdf   ← ape 后端，直接可投稿
  GH7_tree_ggtree.pdf    ← ggtree 后端，ggplot 可继续编辑

着色方案：
  红色  → Flagellate symbiont（鞭毛虫共生体，17条）
  橙色  → Bacteria（细菌，287条）
  紫色  → Fungus outgroup（外类群，8条）
  蓝色  → Other eukaryote（其他真核，198条）
  绿色  → Termite host（白蚁宿主，17条）
  灰色  → Other（未分类，7条）
```

---

## 常见问题

**Q：图里序列名太拥挤看不清？**  
调小 `tip_cex`（ape）或 `tip_label_size`（ggtree），或用 `show_tip_labels = FALSE` 只显示彩色点，另外导出一张带标签的图。

**Q：想去掉 accession 前缀只显示物种名？**  
`classify_tips()` 默认会用 `strip_accession()` 处理 label。如果格式特殊，可自定义 `tip_label_fn` 参数传入你的处理函数。

**Q：内部节点颜色不对？**  
尝试切换 `branch_method = "tip_only"`，只给 tip 着色，内部枝条统一灰色，适合类群边界不清晰的树。
