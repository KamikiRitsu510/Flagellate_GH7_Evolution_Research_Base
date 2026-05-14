# treecolorize

**用 R 给系统发育树自动着色的工具包。**

本包从 GH7 纤维素酶进化研究的 30 余个迭代绘图脚本中提炼而来，核心思路是：给序列名打"关键词标签" → 自动分组 → 自动配色 → 画树。支持 `ape`（静态 PDF/PNG）和 `ggtree`（ggplot2 风格可编辑对象）两种后端。

---

## 安装

```r
# 从本地安装（当前方式）
install.packages("path/to/treecolorize", repos = NULL, type = "source")

# 依赖项
install.packages(c("ape"))
# ggtree 后端（可选）
BiocManager::install(c("ggtree", "treeio"))
install.packages("ggplot2")
```

---

## 最快上手：一行出图

```r
library(ape)
library(treecolorize)

tree <- read.tree("GH7_iqtree.treefile")

# 用内置关键词表自动分类 + 画树，存为 PDF
treecolorize(tree, output_file = "tree_colored.pdf")
```

`treecolorize()` 是一站式封装：内部自动调用分类 → 配色 → 绘图三步。

---

## 分步调用（推荐，灵活度高）

### 第一步：定义分组规则

```r
# 自定义关键词映射表：关键词 → 组名
keyword_map <- list(
  "Reticulitermes"  = "白蚁（低等）",
  "Nasutitermes"    = "白蚁（高等）",
  "Pseudotrichonympha" = "鞭毛虫共生体",
  "Holomastigotoides"  = "鞭毛虫共生体",
  "bacteria"        = "细菌",
  "fungus|Fungi"    = "真菌"      # 支持 | 分隔的多关键词
)
```

> **注意**：关键词大小写不敏感，支持正则表达式。序列名通常是 `accession 物种名` 格式，`classify_tips()` 默认会先去掉 accession 前缀再匹配。

### 第二步：分类

```r
groups <- classify_tips(tree$tip.label, keyword_map)
# 返回命名向量，每个 tip → 组名（未匹配的归入 "Other"）
table(groups)
```

### 第三步：配色

```r
palette <- auto_palette(groups)
# 自动从色盲友好调色板（Wong 2011）取色；超过 8 组时自动扩展
```

### 第四步：画树

```r
# ── ape 后端（推荐，速度快，适合大树）──
plot_tree_ape(tree, groups, palette,
              type = "phylogram",          # "fan" / "unrooted" 同样支持
              branch_method = "majority",  # 枝条着色：按多数子节点颜色
              tip_cex = 0.4,
              title = "GH7 cellulase phylogeny")

# ── 保存为 PDF（自动设置画布尺寸）──
save_tree_ape(tree, groups, "GH7_tree.pdf", palette = palette,
              width = 14, height = 32)

# ── ggtree 后端（返回 ggplot 对象，可继续叠加图层）──
p <- plot_tree_ggtree(tree, groups, palette,
                      layout = "rectangular",
                      tip_label_size = 1.2)
p + geom_treescale()   # 叠加比例尺
ggsave("GH7_tree.pdf", p, width = 14, height = 32)
```

---

## 枝条着色说明

`branch_method` 参数控制内部节点的颜色：

| 选项 | 效果 |
|------|------|
| `"majority"` | 后序遍历，每个内部节点继承直接子节点中占多数的颜色（默认） |
| `"tip_only"` | 仅叶节点着色，内部枝条统一用 `mixed_color` |
| `"uniform"` | 所有枝条同色（`branch_color` 参数指定） |

混合节点（子节点颜色各半）用 `mixed_color`（默认浅灰 `#CCCCCC`）显示。

---

## 常见问题

**Q：序列名里没有物种信息，只有 accession，怎么办？**  
A：先在 R 里替换 `tree$tip.label`，把 accession 映射到物种名再调用本包；或者把物种名直接写进 `keyword_map` 的值而不是键。

**Q：分组太多，颜色区分不开？**  
A：超过 12 个分组时 `auto_palette()` 会自动切换到 `rainbow()`，但视觉效果变差。建议合并次要分组，或手动传入 `palette` 参数指定颜色向量。

**Q：树文件格式不是 Newick，能用吗？**  
A：本包只接受 `ape` 的 `phylo` 对象。`ape::read.nexus()`、`treeio::read.beast()` 等读取后转换一下即可。

---

## 函数速查

| 函数 | 说明 |
|------|------|
| `treecolorize()` | 一站式封装：分类→配色→绘图 |
| `classify_tips()` | 按关键词映射表给 tip 分组 |
| `auto_palette()` | 根据分组数自动生成配色方案 |
| `plot_tree_ape()` | ape 后端绘图（返回不可见，直接出图） |
| `save_tree_ape()` | ape 后端保存 PDF/PNG |
| `plot_tree_ggtree()` | ggtree 后端（返回 ggplot 对象） |
| `branch_colors_majority()` | 计算内部节点颜色（多数投票） |
| `strip_accession()` | 去掉 tip label 里的 accession 前缀 |
| `default_keyword_map()` | 返回内置的常见界/门关键词表 |
