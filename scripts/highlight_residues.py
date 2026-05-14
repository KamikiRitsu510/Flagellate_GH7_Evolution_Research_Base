import py3Dmol

pdb_file = "model_01.pdb"
with open(pdb_file, 'r') as f:
    pdb_data = f.read()

view = py3Dmol.view(width=1200, height=800)
view.addModel(pdb_data, 'pdb')
view.setStyle({'cartoon': {'color': 'gray'}})

residues = [58, 184, 225, 312]
view.addStyle({'resi': residues}, {'stick': {'color': 'red'}, 'sphere': {'color': 'red', 'radius': 0.3}})

view.zoomTo()

html_file = "GH7_positive_sites.html"
with open(html_file, 'w') as f:
    f.write(view._make_html())

print(f"已生成 {html_file}，请用浏览器打开该文件并截图保存。")
