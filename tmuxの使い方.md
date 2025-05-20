# インストール
# 新たなセッションを作る
```bash
tmux new-session -s セッションの名前
```

# 現在のセッションから抜ける

`ctrl+b`を押したあとに`d`

# 作ったセッションを確認
```bash
tmux ls
```

# 作ったセッションを削除
```bash
tmux kill-session -t セッションの名前
```
