#!/usr/bin/env Rscript

# 加载包
library(ape)
library(ggtree)
library(ggplot2)
library(treeio)

# 设置路径（请根据实际位置修改）
old_tree_file <- "11_GH7_tree_consensus.contree"
new_tree_file <- "GH7_final_534.contree"
output_pdf <- "trees_comparison.pdf"

# 读取树
if (!file.exists(old_tree_file)) stop("旧树文件不存在: ", old_tree_file)
if (!file.exists(new_tree_file)) stop("新树文件不存在: ", new_tree_file)

old_tree <- read.tree(old_tree_file)
new_tree <- read.tree(new_tree_file)

cat("旧树叶子数:", length(old_tree$tip.label), "\n")
cat("新树叶子数:", length(new_tree$tip.label), "\n")

# 可选：从旧树中提取新树共有的叶子（如果叶子集不同，建议取消注释）
# common_tips <- intersect(old_tree$tip.label, new_tree$tip.label)
# if (length(common_tips) < length(new_tree$tip.label)) {
#   cat("警告: 新树中有", length(setdiff(new_tree$tip.label, common_tips)), "个叶子不在旧树中，将使用交集\n")
#   old_tree <- keep.tip(old_tree, common_tips)
#   new_tree <- keep.tip(new_tree, common_tips)
# }

# 创建 ggtree 对象（圆形布局）
p1 <- ggtree(old_tree, layout = "circular", branch.length = "branch.length") +
  geom_tiplab(size = 1.5) +
  labs(title = "旧树 (11_GH7_tree_consensus)") +
  theme(plot.title = element_text(hjust = 0.5))

p2 <- ggtree(new_tree, layout = "circular", branch.length = "branch.length") +
  geom_tiplab(size = 1.5) +
  labs(title = "新树 (GH7_final_534)") +
  theme(plot.title = element_text(hjust = 0.5))

# 并排组合
if (!require("patchwork", quietly = TRUE)) install.packages("patchwork")
library(patchwork)
combined <- p1 + p2 + plot_annotation(title = "系统发育树新旧对比")

# 保存 PDF
ggsave(output_pdf, combined, width = 16, height = 8, limitsize = FALSE)
cat("对比图已保存至:", output_pdf, "\n")
