library(ape)

# 读入树
tree <- read.tree("results/GH7_final_534.contree")

# 鞭毛虫关键词
flag_patterns <- c("AAY833", "BAB645", "BAC075")
flag_labels <- tree$tip.label[grepl(paste(flag_patterns, collapse="|"), tree$tip.label)]

# 水蚤外群
daphnia_labels <- tree$tip.label[grepl("Daphnia|pulex", tree$tip.label, ignore.case=TRUE)]

# 真菌外群（指定登录号）
fungi_acc <- c("XP_006965674.1", "XP_003054277.1", "AGT15838.1")
fungi_labels <- fungi_acc[fungi_acc %in% tree$tip.label]

# 合并并去重
keep <- unique(c(flag_labels, daphnia_labels, fungi_labels))

# 剪枝
subtree <- drop.tip(tree, tree$tip.label[!(tree$tip.label %in% keep)])

# 保存
write.tree(subtree, file = "results/subtree_flagellate.nwk")

cat("原始序列数:", length(tree$tip.label), "\n")
cat("保留序列数:", length(subtree$tip.label), "\n")
cat("子树已保存: results/subtree_flagellate.nwk\n\n保留的序列:\n")
print(keep)
