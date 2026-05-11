#!/usr/bin/env Rscript

# 加载包
library(ape)
library(phytools)

# 文件路径
old_file <- "11_GH7_tree_consensus.contree"
new_file <- "GH7_final_534.contree"
out_pdf <- "trees_comparison_ape.pdf"

# 检查文件是否存在
if (!file.exists(old_file)) stop("旧树文件不存在: ", old_file)
if (!file.exists(new_file)) stop("新树文件不存在: ", new_file)

# 读树
old_tree <- read.tree(old_file)
new_tree <- read.tree(new_file)

cat("旧树叶子数:", length(old_tree$tip.label), "\n")
cat("新树叶子数:", length(new_tree$tip.label), "\n")

# 可选：如果叶子数差异太大，为了避免重叠，可以只保留共有叶子
common_tips <- intersect(old_tree$tip.label, new_tree$tip.label)
if (length(common_tips) < length(new_tree$tip.label)) {
  cat("警告：新树中有", length(setdiff(new_tree$tip.label, common_tips)), "个叶子不在旧树中。\n")
  cat("将使用交集（", length(common_tips), "个叶子）进行对比。\n")
  old_tree <- keep.tip(old_tree, common_tips)
  new_tree <- keep.tip(new_tree, common_tips)
}

# 打开 PDF
pdf(out_pdf, width = 12, height = 8)

# 并排绘制
par(mfrow = c(1, 2), mar = c(1, 1, 3, 1))

# 旧树（圆形布局）
plotTree(old_tree, type = "fan", fsize = 0.6, lwd = 1, 
         ftype = "i", offset = 0.5)
title("旧树 (11_GH7_tree_consensus)", cex.main = 1.2)

# 新树（圆形布局）
plotTree(new_tree, type = "fan", fsize = 0.6, lwd = 1,
         ftype = "i", offset = 0.5)
title("新树 (GH7_final_534)", cex.main = 1.2)

dev.off()
cat("对比图已保存至:", out_pdf, "\n")
