# Doxygen MCP Server

Doxygenが生成したFrontISTRのXMLドキュメントをAIアシスタント（GitHub Copilot、Claude等）が検索・分析できるようにします。

## 機能

- 関数・サブルーチン・変数のキーワード検索
- モジュール・クラス・構造体の検索
- 関数呼び出し関係の追跡（caller/callee解析）
- ソースコード位置の特定

## セットアップ

### 要件

- Python 3.7以上
- Doxygen
- fastmcp

### インストール

Python仮想環境を作成し、必要なパッケージをインストールします：

```bash
# 仮想環境の作成
python -m venv ai-tools/.venv

# 仮想環境の有効化
# Linux/macOS:
source ai-tools/.venv/bin/activate
# Windows (PowerShell):
ai-tools\.venv\Scripts\Activate.ps1
# Windows (cmd):
ai-tools\.venv\Scripts\activate.bat

# パッケージのインストール
pip install fastmcp
```

### Doxygenドキュメントの生成

CMakeビルドディレクトリで `-DWITH_DOC=ON` を有効にしてドキュメントを生成します：

```bash
cd build
cmake -DWITH_DOC=ON ..
make doc
```

これにより `build/doc/xml/` にXMLドキュメントが生成されます。

### MCP設定

#### VS Code（GitHub CopilotまたはClaude拡張機能）の場合

`.vscode/mcp.json` に以下の設定を記述し、VS Codeを再起動してMCPサーバーを起動します。

**Linux/macOS:**
```json
{
  "servers": {
    "doxygen": {
      "type": "stdio",
      "command": "${workspaceFolder}/ai-tools/.venv/bin/python",
      "args": [
        "${workspaceFolder}/ai-tools/doxygen_mcp.py",
        "--xml", "${workspaceFolder}/build/doc/xml",
        "--html", "${workspaceFolder}/build/doc/html"
      ]
    }
  }
}
```

**Windows:**
```json
{
  "servers": {
    "doxygen": {
      "type": "stdio",
      "command": "${workspaceFolder}/ai-tools/.venv/Scripts/python.exe",
      "args": [
        "${workspaceFolder}/ai-tools/doxygen_mcp.py",
        "--xml", "${workspaceFolder}/build/doc/xml",
        "--html", "${workspaceFolder}/build/doc/html"
      ]
    }
  }
}
```

#### その他のAIアシスタント

各ツールのMCP設定に従って、以下のコマンドを登録してください：

**Linux/macOS:**
```bash
ai-tools/.venv/bin/python ai-tools/doxygen_mcp.py --xml build/doc/xml --html build/doc/html
```

**Windows:**
```cmd
ai-tools\.venv\Scripts\python.exe ai-tools\doxygen_mcp.py --xml build\doc\xml --html build\doc\html
```

## 使用方法

GitHub Copilot ChatまたはClaude等のAIアシスタントに自然言語で質問します。MCPサーバーは自動的に呼び出されます。

**クエリ例**：
- "fstr_newton関数の呼び出し先を列挙してください"
- "動的解析でダンピング行列を扱う関数を検索してください"
- "m_stepモジュールの関数一覧を取得してください"

