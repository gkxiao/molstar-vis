# Mol* 可视化交付物生成器

将 PDB 格式的受体结构文件与对接后的 SDF 格式配体文件，打包成一个可在浏览器中直接打开的 Mol* 可视化交付物（`molstar_viewer.zip`）。

## 项目结构

```
├── 4ivc_receptor.pdb          # PDB 格式受体结构文件
├── test.sdf                   # 对接好的配体文件（V2000 SDF 格式）
├── prep_data.py               # 准备可视化交付物的 Python 脚本
├── molstar_viewer.zip         # 可视化交付物模板（解压后为 molstar_viewer 目录）
└── README.md
```

解压 `molstar_viewer.zip` 后得到模板目录：

```
molstar_viewer/
├── data.js                    # 数据文件（由脚本生成）
├── index.html                 # 可视化页面入口
└── lib/
    ├── molstar-bundle.js      # Mol* 核心库
    ├── molstar.css            # Mol* 样式
    └── MOLSTAR_LICENSE.txt    # Mol* 许可证
```

用浏览器打开 `index.html` 即可在 Mol* 中可视化结果。

## 环境要求

- Python 3
- 依赖以下标准库模块（无需额外安装）：

  ```
  argparse, glob, json, re, shutil, subprocess, sys
  collections, pathlib
  ```

  以及脚本内使用的其他模块（视实际代码而定）。如脚本无第三方依赖，可直接用系统 Python 运行。

## 使用方法

1. 克隆本仓库：

   ```bash
   git clone https://github.com/gkxiao/molstar-vis.git
   cd molstar-vis
   ```

2. 确保以下内容处于**同一目录**下：
   - `prep_data.py`
   - PDB 格式的受体文件（如 `4ivc_receptor.pdb`）
   - 对接好的配体 SDF 文件（如 `test.sdf`，V2000 版本）
   - 解压后的 `molstar_viewer` 目录

3. 解压模板（若尚未解压）：

   ```bash
   unzip molstar_viewer.zip
   ```

4. 运行脚本生成交付物：

   ```bash
   python prep_data.py
   ```

5. 脚本运行结束后，当前目录下会生成 `molstar_viewer.zip`。解压后打开 `index.html` 即可在 Mol* 中查看可视化结果。

## 输出

- `molstar_viewer.zip`：包含受体结构、配体对接结果与 Mol* 可视化页面的完整交付包。
- 解压后目录结构如上，浏览器打开 `index.html` 即可交互查看。

## 说明

- 配体 SDF 文件需为 **V2000** 格式，以确保 Mol* 正确解析。
- 受体 PDB 文件与配体 SDF 文件应处于同一坐标系（即对接后的坐标），以保证可视化时二者正确叠合。
- `data.js` 由 `prep_data.py` 自动生成，请勿手动修改。

## 许可证

- 本项目代码遵循仓库内声明的许可证。
- Mol* 相关文件（`lib/` 下）遵循其自身许可证，详见 `MOLSTAR_LICENSE.txt`。

## 参考

- [Mol* 官方网站](https://molstar.org/)
- [Mol* GitHub 仓库](https://github.com/molstar/molstar)
