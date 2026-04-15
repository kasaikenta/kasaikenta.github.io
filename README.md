# kasaikenta.github.io

This GitHub Pages repository is now a redirect-only compatibility endpoint.
The canonical homepage is:

<https://kasai.ict.eng.isct.ac.jp/>

All tracked HTML pages in this repository should stay as lightweight redirect
pages that preserve the path, query string, and fragment identifier. The
canonical content source is the local campus-site tree:

`/Users/kasaikenta/Projects/kasaikenta_webpage`

Do not restore full page contents in this repository unless intentionally
rolling back the campus-site migration. Static assets may remain here for a
while to avoid breaking old direct links.

## Google検索キーワードの確認方法（GA4）

GA4/Search Console の本体確認対象は `https://kasai.ict.eng.isct.ac.jp/`
です。GitHub Pages 側はリダイレクト専用なので、検索・アクセス解析の
正規URLも学内サイト側に寄せます。

Google検索キーワード（どの語句で流入したか）を見るには、GA4 と Search Console の連携が必要です。

1. Search Console でサイトプロパティを作成・確認する
2. GA4 管理画面で `管理 > プロダクトのリンク > Search Console のリンク` から連携する
3. GA4 の `レポート > ライブラリ` で Search Console コレクションを公開する
4. `レポート > Search Console > クエリ` で検索語句を確認する

注意:
- データ反映には通常 24〜48 時間程度かかります。
- ホームページ（`/`）に限定した検索語句を厳密に見る場合は、Search Console 側の `検索パフォーマンス` でページを `/` にフィルタして確認するのが確実です。
