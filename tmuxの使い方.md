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
# 作ったセッションに入る
```bash
tmux attach-session -t <セッション名または番号>
```

# 作ったセッションを削除
```bash
tmux kill-session -t セッションの名前
```
