// docs/.vitepress/theme/redirects.ts: moved pages of the docs, old site path -> new site path (relative to the site
// base, without .html), the scheme of ADAM's docs theme. GitHub Pages has no server redirects, so an old URL renders the
// 404 page and index.ts replaces it, client side, by the new location. Keep the entries: external links and bookmarks
// still point at the old locations.

export const pageRedirects: Record<string, string> = {
  // 2026-11-15: the guide's short Verification page merged into the Verification section's landing page
  '/guide/verification': '/verification/',
}
