# GitHub 零基础入门指南（中文版）

> 适用对象：完全没有使用过 Git / GitHub 的同学。
> 目标：把本仓库（`las-crc-prediction`）上传到自己的 GitHub 账号，并学会后续更新。
> ⚠️ **最重要的安全警示（先读这一段）**
> - **绝不上传 TCGA 原始数据或任何患者级数据**（`data/raw/`、`data/processed/` 下的 rds / csv 文件）。
>   TCGA 数据受 GDC 数据使用协议约束，不允许公开再分发；患者数据涉及隐私与伦理。
> - 本仓库的 `.gitignore` 已经配置好，会自动忽略这些目录和文件，但你在执行
>   `git add` 前仍应检查 `git status`，确认没有数据文件被列进去。
> - 一旦误传，文件会留在 Git 历史中，仅删除当前文件**不够**，需要重写历史；
>   预防远比补救重要。

---

## 第 1 步：注册 GitHub 账号

1. 打开 <https://github.com>，点击右上角 **Sign up**。
2. 填写邮箱、密码、用户名（用户名将出现在仓库网址中，建议用真实姓名拼音或学术化名，如 `xingxing-tao`）。
3. 完成邮箱验证。
4. （可选）在 <https://education.github.com> 用学校邮箱申请 GitHub Education 学生包，可获 Copilot 免费使用权等福利。

## 第 2 步：安装 Git

- **Windows**：到 <https://git-scm.com/download/win> 下载安装包，一路默认下一步即可。安装后右键任意文件夹应能看到 “Git Bash Here”。
- **macOS**：打开“终端”（Terminal），输入 `git --version`；若未安装，系统会提示安装 Xcode 命令行工具，按提示操作。
- **Linux**：`sudo apt install git`（Debian/Ubuntu）或 `sudo dnf install git`（Fedora）。

安装完成后验证：

```bash
git --version    # 显示版本号即安装成功，例如 git version 2.43.0
```

## 第 3 步：首次配置 Git（全局，只需一次）

把下面的名字和邮箱换成你自己的（邮箱建议与 GitHub 注册邮箱一致）：

```bash
git config --global user.name "Xingxing Tao"
git config --global user.email "your_email@example.com"
```

检查配置：

```bash
git config --global --list
```

## 第 4 步：关于本仓库的 git 初始化

本仓库目录里**已经执行过** `git init` 并完成了第一次提交（commit），
你不需要再初始化。进入目录即可查看：

```bash
cd las-crc-prediction
git log --oneline     # 应能看到一条 "Initial commit: ..." 的记录
git status            # 应显示 working tree clean
```

## 第 5 步：在 GitHub 网页端创建远程仓库

1. 登录 GitHub，点击右上角 **+** → **New repository**。
2. Repository name 填 `las-crc-prediction`（与本目录同名，方便对应）。
3. 建议选 **Public**（作为申请材料展示用；如暂不想公开可选 Private，之后再改）。
4. **不要勾选** “Add a README file”、“Add .gitignore”、“Choose a license”
   ——本仓库已包含这些文件，勾选会导致首次推送冲突。
5. 点击 **Create repository**。

创建后页面会显示仓库地址，形如：

```
https://github.com/<你的用户名>/las-crc-prediction.git
```

## 第 6 步：把本地仓库关联到远程并推送

在本仓库目录下执行（把 `<你的用户名>` 替换掉）：

```bash
git remote add origin https://github.com/<你的用户名>/las-crc-prediction.git
git branch -M main          # 把主分支命名为 main（GitHub 现在的默认约定）
git push -u origin main     # 第一次推送；-u 让以后直接 git push 即可
```

首次推送会要求登录认证：

- 推荐方式：安装 [GitHub CLI](https://cli.github.com/) 后运行 `gh auth login`，按提示用浏览器授权；
- 或者按提示输入用户名 + **Personal Access Token**（注意：GitHub 已不接受账号密码，
  Token 在 GitHub 网页 Settings → Developer settings → Personal access tokens → Tokens (classic) 生成，勾选 `repo` 权限）。

推送成功后，刷新 GitHub 仓库网页，就能看到全部文件了。

## 第 7 步：日常更新流程（改了代码之后）

```bash
git status                          # 1. 查看哪些文件变了（确认没有数据文件！）
git add README.md R/03_las_ssgsea.R # 2. 添加指定文件；或 git add -A 添加全部
git commit -m "Add construct validity check for LAS"   # 3. 提交，写清改动说明
git push                            # 4. 推送到 GitHub
```

提交信息（commit message）建议用英文、动词开头、一句话说清做了什么。

## 备选方案：网页端直接上传（不推荐，仅应急）

GitHub 仓库页面点 **Add file → Upload files**，可以拖拽上传单个文件。
缺点：不能上传空目录、容易遗漏 `.gitignore`、没有版本历史记录。
仅在完全无法用命令行时应急使用，正式使用请学习上面的命令行流程。

## 常见问题

| 问题 | 解决办法 |
|------|----------|
| `push` 时提示 `rejected: fetch first` | 远程有本地没有的提交，先 `git pull --rebase` 再 push |
| 误把数据文件 `git add` 了但还没 commit | `git restore --staged <文件>` 撤回到未暂存状态 |
| 误 commit 了数据文件但还没 push | `git rm --cached <文件>` 后重新 commit |
| 已 push 了不该传的数据 | 立即在 GitHub 把仓库设为 Private，并参考 GitHub 官方文档 “Removing sensitive data from a repository” 重写历史 |
| 大文件（>100 MB）传不上 | GitHub 单文件上限 100 MB；TCGA 数据本就不该上传，确认它在 `.gitignore` 中 |

## 上传前最后检查清单

- [ ] `git status` 中没有 `data/raw/`、`data/processed/`、`results/` 下的任何文件
- [ ] 没有任何 rds / 数据 csv 进入版本库
- [ ] README.md 中的占位信息（如联系方式）已替换为真实内容
- [ ] 仓库可见性（Public/Private）符合你的预期
