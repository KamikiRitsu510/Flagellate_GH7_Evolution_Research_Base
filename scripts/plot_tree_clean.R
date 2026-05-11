#!/usr/bin/env Rscript

# ======================================================
# 绘制最终系统发育树（矩形，按物种着色，显示物种名）
# ======================================================

# 安装缺失的包（如果已安装则跳过）
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos = "https://cran.rstudio.com")
if (!requireNamespace("ggtree", quietly = TRUE))
  BiocManager::install("ggtree", update = FALSE, ask = FALSE)
if (!requireNamespace("ggplot2", quietly = TRUE))
  install.packages("ggplot2", repos = "https://cran.rstudio.com")
if (!requireNamespace("treeio", quietly = TRUE))
  BiocManager::install("treeio", update = FALSE, ask = FALSE)

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggtree)
  library(treeio)
})

# ------------- 用户配置 -------------
tree_file <- "GH7_final_534.contree"   # 请确保此文件存在
output_pdf <- "GH7_tree_clean.pdf"

# ------------- 读取树 -------------
if (!file.exists(tree_file)) stop("树文件不存在: ", tree_file)
tree <- read.tree(tree_file)

# ------------- 从 tip 标签中提取物种名和分类 -------------
# 假设你的 fasta 标题格式为: "登录号 [空格] 物种名"
# 例如: "AAY83390.3 Pseudotrichonympha grassii"
# 如果格式不同，请修改下面的正则表达式
tip_labels <- tree$tip.label
species <- sub("^[A-Za-z0-9_.]+\\s+", "", tip_labels)   # 去掉登录号，保留物种名
# 如果物种名中包含多余信息（比如 strain），可以进一步清理，但先保留完整

# 定义分类规则（基于物种名中的关键词）
classify <- function(sp) {
  sp_lower <- tolower(sp)
  if (grepl("flagellate|pseudotrichonympha|holomastigotoides|trichonympha|spirotrichonympha", sp_lower))
    return("Flagellate")
  if (grepl("bacteria|bacillus|clostridium|proteobacteria|actinobacteria", sp_lower))
    return("Bacteria")
  if (grepl("fungus|fungi|ascomycota|basidiomycota|chaetomium|phytophthora|pyrenophora", sp_lower))
    return("Fungi")
  if (grepl("daphnia|water flea|crustacea", sp_lower))
    return("Crustacea")
  return("Other")
}
group <- sapply(species, classify)

# 颜色方案
col_pal <- c("Flagellate" = "#E41A1C", "Bacteria" = "#377EB8", 
             "Fungi" = "#4DAF4A", "Crustacea" = "#984EA3", 
             "Other" = "#FF7F00")

# 构建数据框用于着色
df <- data.frame(label = tip_labels, group = group, species = species)

# ------------- 绘制树（矩形布局，非螺旋）-------------
p <- ggtree(tree, layout = "rectangular", branch.length = "branch.length") %<+% df +
  geom_tippoint(aes(color = group), size = 0.8) +
  geom_tiplab(aes(label = species), size = 1.5, offset = 0.5, align = TRUE, linesize = 0.2) +
  scale_color_manual(values = col_pal, name = "Group") +
  theme(legend.position = "right",
        plot.margin = unit(c(0.5, 2, 0.5, 0.5), "cm")) +
  labs(title = "GH7 Phylogenetic Tree (534 sequences)")

# 保存 PDF（宽高可调）
ggsave(output_pdf, p, width = 14, height = 20, limitsize = FALSE)
cat("✅ 树图已保存至:", output_pdf, "\n")
cat("总序列数:", length(tip_labels), "\n")
cat("分类统计:\n")
print(table(group))
