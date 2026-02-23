# kasaikenta.github.io

## Google検索キーワードの確認方法（GA4）

このリポジトリの `index.html` には GA4 タグ（`G-YY8X0H44J2`）を設定済みです。

Google検索キーワード（どの語句で流入したか）を見るには、GA4 と Search Console の連携が必要です。

1. Search Console でサイトプロパティを作成・確認する
2. GA4 管理画面で `管理 > プロダクトのリンク > Search Console のリンク` から連携する
3. GA4 の `レポート > ライブラリ` で Search Console コレクションを公開する
4. `レポート > Search Console > クエリ` で検索語句を確認する

注意:
- データ反映には通常 24〜48 時間程度かかります。
- ホームページ（`/`）に限定した検索語句を厳密に見る場合は、Search Console 側の `検索パフォーマンス` でページを `/` にフィルタして確認するのが確実です。
