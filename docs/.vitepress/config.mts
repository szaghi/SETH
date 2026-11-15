import { withMermaid } from 'vitepress-plugin-mermaid'

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
      { text: 'Home',   link: '/' },
      { text: 'Guide',  link: '/guide/' },
      { text: 'Verification', link: '/verification/' },
      { text: 'GitHub', link: 'https://github.com/szaghi/SETH' },
    ],

    // Add your guide pages here.
    sidebar: {
      '/verification/': [
        { text: 'Verification', items: [
          { text: 'Catalog',       link: '/verification/' },
          { text: 'Test programs', link: '/verification/unit-tests' },
        ] },
        { text: 'Accuracy', collapsed: false, items: [
          { text: 'Overview',                    link: '/verification/accuracy/' },
          { text: 'A1 Order of accuracy',        link: '/verification/accuracy/order-of-accuracy' },
          { text: 'A2 Critical points',          link: '/verification/accuracy/critical-points' },
          { text: 'A3 Whole-domain convergence', link: '/verification/accuracy/convergence' },
          { text: 'A4 Euler order of accuracy',  link: '/verification/accuracy/euler-order' },
        ] },
        { text: 'Shock capturing', collapsed: false, items: [
          { text: 'Overview',                  link: '/verification/shocks/' },
          { text: 'S1 Sod',                    link: '/verification/shocks/sod' },
          { text: 'S2 Lax',                    link: '/verification/shocks/lax' },
          { text: 'S3 123 problem',            link: '/verification/shocks/test123' },
          { text: 'S4 Shu–Osher',              link: '/verification/shocks/shu-osher' },
          { text: 'S5 Woodward–Colella blast', link: '/verification/shocks/blast' },
          { text: 'S6 Burgers',                link: '/verification/shocks/burgers' },
          { text: 'S7 Composite profile',      link: '/verification/shocks/advection' },
        ] },
        { text: 'Spectral', collapsed: false, items: [
          { text: 'Overview',                  link: '/verification/spectral/' },
          { text: 'F1 Linear backgrounds',     link: '/verification/spectral/linear-spectra' },
          { text: 'F2 Dispersion relation',    link: '/verification/spectral/adr' },
        ] },
        { text: 'Contracts', collapsed: false, items: [
          { text: 'Overview',                  link: '/verification/contracts/' },
          { text: 'C1 Equivalence with ADAM',  link: '/verification/contracts/adam-equivalence' },
          { text: 'C2 Properties',             link: '/verification/contracts/properties' },
          { text: 'C3 Coefficient tables',     link: '/verification/contracts/tables' },
          { text: 'C4 Golden regression',      link: '/verification/contracts/golden' },
          { text: 'C5 Device',                 link: '/verification/contracts/device' },
          { text: 'C6 Mutation testing',       link: '/verification/contracts/mutation' },
          { text: 'C7 Positivity limiter',     link: '/verification/contracts/positivity' },
        ] },
      ],
      '/guide/': [
        {
          text: 'Introduction',
          items: [
            { text: 'About',        link: '/guide/' },
            { text: 'Installation', link: '/guide/installation' },
            { text: 'Quick Start',  link: '/guide/quickstart' },
          ],
        },
        {
          text: 'Reference',
          items: [
            { text: 'Conventions',  link: '/guide/conventions' },
            { text: 'Schemes',      link: '/guide/schemes' },
            { text: 'API',          link: '/guide/api' },
          ],
        },
        {
          text: 'Quality',
          items: [
            { text: 'Verification', link: '/guide/verification' },
            { text: 'Performance',  link: '/guide/performance' },
            { text: 'Sources',      link: '/guide/sources' },
          ],
        },
        {
          text: 'Project',
          items: [
            { text: 'Contributing', link: '/guide/contributing' },
            { text: 'Changelog',    link: '/guide/changelog' },
          ],
        },
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
