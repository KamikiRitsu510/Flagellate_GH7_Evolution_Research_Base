from ete3 import Tree

# 加载树文件
t = Tree("ASR_iqtree.treefile")

# 17条鞭毛虫序列的登录号
target_taxa = [
    "BAB64565.3", "BAC07551.1", "BAC07552.1", "AAY83390.3", "BAB64553.1",
    "BAB64554.1", "BAB64555.1", "BAB64556.1", "BAB64557.1", "BAB64558.1",
    "BAB64559.1", "BAB64560.1", "BAB64561.1", "BAB64562.1", "BAB64563.2",
    "BAB64564.2", "AAY83391.1"
]

# 计算并打印MRCA节点
mrca = t.get_common_ancestor(target_taxa)
print(f"这些序列的最近共同祖先(MRCA)节点是: {mrca.name}")
