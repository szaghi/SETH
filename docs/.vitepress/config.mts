import { existsSync, readFileSync } from 'node:fs'
import { resolve } from 'node:path'
import { withMermaid } from 'vitepress-plugin-mermaid'

// formal's sidebar of the API pages (docs/api/_sidebar.json), when the API has been generated
const apiSidebarFile = resolve(__dirname, '../api/_sidebar.json')
const apiSidebar = existsSync(apiSidebarFile) ? JSON.parse(readFileSync(apiSidebarFile, 'utf-8')) : []

export default withMermaid({
  title: 'SETH',
  description: 'Smoothness-weighted ENO, TENO and Hybrids: high-order non-oscillatory reconstruction and interpolation for HPC',
  base: '/SETH/',
  srcExclude: ['ford.md', '**/_tables/**'],

  // the identity (Ogi in Neon Moegi): favicons by colour scheme, touch icon, theme colour
  head: [
    ['link', { rel: 'icon', type: 'image/svg+xml', href: '/SETH/favicon.svg', media: '(prefers-color-scheme: light)' }],
    ['link', { rel: 'icon', type: 'image/svg+xml', href: '/SETH/favicon-dark.svg', media: '(prefers-color-scheme: dark)' }],
    ['link', { rel: 'icon', type: 'image/png', sizes: '32x32', href: '/SETH/favicon-32.png' }],
    ['link', { rel: 'icon', type: 'image/png', sizes: '16x16', href: '/SETH/favicon-16.png' }],
    ['link', { rel: 'apple-touch-icon', sizes: '180x180', href: '/SETH/apple-touch-icon.png' }],
    ['meta', { name: 'theme-color', content: '#f3ede1', media: '(prefers-color-scheme: light)' }],
    ['meta', { name: 'theme-color', content: '#0e1116', media: '(prefers-color-scheme: dark)' }],
    // the theme-aware videos of the verification pages: <video class="light-only"> / <video class="dark-only">
    ['style', {}, '.dark-only{display:none}.dark .dark-only{display:block}.dark .light-only{display:none}'],
  ],

  markdown: {
    math: true,
    languages: ['fortran-free-form', 'fortran-fixed-form'],
    languageAlias: {
      fortran: 'fortran-free-form',
      f90: 'fortran-free-form',
      f03: 'fortran-free-form',
      f08: 'fortran-free-form',
    },
  },

  themeConfig: {
    logo: { light: '/logo-nav.svg', dark: '/logo-nav-dark.svg', alt: 'SETH' },
    siteTitle: 'SETH',
    nav: [
      { text: 'Start',        link: '/guide/' },
      { text: 'Tutorial',     link: '/tutorial/' },
      { text: 'Cookbook',     link: '/cookbook/' },
      { text: 'Reference',    link: '/reference/' },
      { text: 'Verification', link: '/verification/' },
      { text: 'API',          link: '/api/' },
      { text: 'GitHub',       link: 'https://github.com/szaghi/SETH' },
    ],

    // One sidebar for the whole site, in reading order: it is also the chain of the Previous / Next links at the foot
    // of every page, so reading straight through goes start -> tutorial -> recipes -> reference -> verification.
    // A collapsed group opens by itself on its own pages. scripts/check_docs.py fails when a page is missing here.
    sidebar: {
      // the module pages written by formal (docs:build runs it first), with its own sidebar
      '/api/': [
        { text: 'Back to the documentation', link: '/guide/' },
        { text: 'API', items: [{ text: 'Overview', link: '/api/' }] },
        ...apiSidebar,
      ],
      '/': [
        { text: 'Start here', items: [
          { text: 'About',        link: '/guide/' },
          { text: 'Installation', link: '/guide/installation' },
          { text: 'Quick start',  link: '/guide/quickstart' },
          { text: 'Conventions',  link: '/guide/conventions' },
        ] },
        { text: 'Tutorial', collapsed: false, items: [
          { text: 'All chapters',            link: '/tutorial/' },
          { text: '1. One face by hand',     link: '/tutorial/01-face' },
          { text: '2. A line of faces',      link: '/tutorial/02-line' },
          { text: '3. A conservation law',   link: '/tutorial/03-advection' },
          { text: '4. The weights',          link: '/tutorial/04-weights' },
          { text: '5. TENO and its cutoff',  link: '/tutorial/05-teno' },
          { text: '6. A nonlinear flux',     link: '/tutorial/06-burgers' },
          { text: '7. A system',             link: '/tutorial/07-euler' },
          { text: '8. Positivity',           link: '/tutorial/08-positivity' },
          { text: '9. A hybrid',             link: '/tutorial/09-hybrid' },
          { text: '10. Non-polynomial TENO', link: '/tutorial/10-nonpoly' },
          { text: '11. On the GPU',          link: '/tutorial/11-gpu' },
        ] },
        { text: 'Cookbook', collapsed: false, items: [
          { text: 'All recipes',                   link: '/cookbook/' },
          { text: 'Choose a scheme',               link: '/cookbook/choose' },
          { text: 'Point values or cell averages', link: '/cookbook/point-values' },
          { text: 'Flag the shocks',               link: '/cookbook/shock-flags' },
          { text: 'Bound a reconstruction',        link: '/cookbook/mp-bound' },
          { text: 'A central branch for Euler',    link: '/cookbook/central' },
          { text: 'Single precision',              link: '/cookbook/precision' },
          { text: 'Time a kernel on a GPU',        link: '/cookbook/performance' },
        ] },
        { text: 'Reference', collapsed: false, items: [
          { text: 'Feature map',         link: '/reference/' },
          { text: 'WENO-JS',             link: '/reference/weno-js' },
          { text: 'WENO-Z',              link: '/reference/weno-z' },
          { text: 'TENO',                link: '/reference/teno' },
          { text: 'Centred WENO',        link: '/reference/weno-c' },
          { text: 'WENO-6 point interpolation', link: '/reference/weno6-point' },
          { text: 'Non-polynomial TENO', link: '/reference/nonpolynomial' },
          { text: 'Linear schemes',      link: '/reference/linear' },
          { text: 'Sensors',             link: '/reference/sensors' },
          { text: 'Central fluxes',      link: '/reference/central' },
          { text: 'Limiters',            link: '/reference/limiters' },
          { text: 'Kinds and constants', link: '/reference/kinds' },
          { text: 'Schemes: the theory', link: '/guide/schemes' },
          { text: 'Precision',           link: '/guide/precision' },
          { text: 'Performance',         link: '/guide/performance' },
        ] },
        { text: 'Verification', collapsed: true, items: [
          { text: 'Catalog',       link: '/verification/' },
          { text: 'Test programs', link: '/verification/unit-tests' },
        ] },
        { text: 'Verification: Accuracy', collapsed: true, items: [
          { text: 'Overview',                    link: '/verification/accuracy/' },
          { text: 'A1 Order of accuracy',        link: '/verification/accuracy/order-of-accuracy' },
          { text: 'A2 Critical points',          link: '/verification/accuracy/critical-points' },
          { text: 'A3 Whole-domain convergence', link: '/verification/accuracy/convergence' },
          { text: 'A4 Euler order of accuracy',  link: '/verification/accuracy/euler-order' },
          { text: 'A5 Centred WENO',             link: '/verification/accuracy/centred-weno' },
        ] },
        { text: 'Verification: Shock capturing', collapsed: true, items: [
          { text: 'Overview',                  link: '/verification/shocks/' },
          { text: 'S1 Sod',                    link: '/verification/shocks/sod' },
          { text: 'S2 Lax',                    link: '/verification/shocks/lax' },
          { text: 'S3 123 problem',            link: '/verification/shocks/test123' },
          { text: 'S4 Shu–Osher',              link: '/verification/shocks/shu-osher' },
          { text: 'S5 Woodward–Colella blast', link: '/verification/shocks/blast' },
          { text: 'S6 Burgers',                link: '/verification/shocks/burgers' },
          { text: 'S7 Composite profile',      link: '/verification/shocks/advection' },
        ] },
        { text: 'Verification: Spectral', collapsed: true, items: [
          { text: 'Overview',                  link: '/verification/spectral/' },
          { text: 'F1 Linear backgrounds',     link: '/verification/spectral/linear-spectra' },
          { text: 'F2 Dispersion relation',    link: '/verification/spectral/adr' },
        ] },
        { text: 'Verification: Hybrids', collapsed: true, items: [
          { text: 'Overview',                  link: '/verification/hybrids/' },
          { text: 'H1 Sensors',                link: '/verification/hybrids/sensors' },
          { text: 'H2 Central fluxes',         link: '/verification/hybrids/central' },
          { text: 'H3 Hybrids in a solver',    link: '/verification/hybrids/harness' },
        ] },
        { text: 'Verification: Non-polynomial', collapsed: true, items: [
          { text: 'Overview',                  link: '/verification/nonpolynomial/' },
          { text: 'NP1 Flags',                 link: '/verification/nonpolynomial/flags' },
          { text: 'NP2 Detectors',             link: '/verification/nonpolynomial/detectors' },
          { text: 'NP3 THINC override',        link: '/verification/nonpolynomial/thinc' },
          { text: 'NP4 DR override',           link: '/verification/nonpolynomial/dr' },
          { text: 'NP5 Stencil targeting',     link: '/verification/nonpolynomial/targeting' },
          { text: 'NP6 Gaussian-RBF fluxes',   link: '/verification/nonpolynomial/rbf' },
          { text: 'NP7 Ablation',              link: '/verification/nonpolynomial/ablation' },
          { text: 'NP8 Length units',          link: '/verification/nonpolynomial/units' },
          { text: 'NP9 Two dimensions',        link: '/verification/nonpolynomial/two-d' },
        ] },
        { text: 'Verification: Contracts', collapsed: true, items: [
          { text: 'Overview',                  link: '/verification/contracts/' },
          { text: 'C1 Equivalence with ADAM',  link: '/verification/contracts/adam-equivalence' },
          { text: 'C2 Properties',             link: '/verification/contracts/properties' },
          { text: 'C3 Coefficient tables',     link: '/verification/contracts/tables' },
          { text: 'C4 Golden regression',      link: '/verification/contracts/golden' },
          { text: 'C5 Device',                 link: '/verification/contracts/device' },
          { text: 'C6 Mutation testing',       link: '/verification/contracts/mutation' },
          { text: 'C7 Positivity limiter',     link: '/verification/contracts/positivity' },
        ] },
        { text: 'Project', collapsed: true, items: [
          { text: 'Tooling',      link: '/guide/tooling' },
          { text: 'Sources',      link: '/guide/sources' },
          { text: 'Contributing', link: '/guide/contributing' },
          { text: 'Changelog',    link: '/guide/changelog' },
        ] },
      ],
    },

    socialLinks: [
      { icon: 'github', link: 'https://github.com/szaghi/SETH' },
    ],

    search: {
      provider: 'local',
    },

    footer: {
      message: 'Released under your choice of the GPL v3, BSD-2-Clause, BSD-3-Clause or MIT license.',
      copyright: 'Copyright © 2026 Stefano Zaghi',
    },
  },

  mermaid: {},

  vite: {
    // Build with an explicit modern JS target so the docs compile regardless of
    // which mermaid/vitepress/esbuild versions npm resolves. Vite's default
    // es2020 target forces esbuild to down-level modern syntax (e.g. the
    // destructuring mermaid 11.16+ emits), which it refuses to do and the build
    // dies. es2022 needs no lowering and is within VitePress's browser floor.
    build: {
      target: 'es2022',
    },
    // The dev server pre-bundles the dependencies with its own esbuild target (es2020), ignoring build.target:
    // without this, `vitepress dev` dies on the destructuring of mark.js and @vueuse. Bundling mermaid whole also
    // covers the CJS dependencies vitepress-plugin-mermaid does not list (fastdom since mermaid 11.17), which leave
    // a white dev page otherwise.
    optimizeDeps: {
      include: ['mermaid'],
      esbuildOptions: {
        target: 'es2022',
      },
    },
  },
})
