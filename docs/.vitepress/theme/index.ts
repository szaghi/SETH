// docs/.vitepress/theme/index.ts
// The SETH identity (custom.css: the Sumi family tokens in the Neon Moegi ink, Ogi mark) plus the dark figures of the ADAM docs
// harness (copied from ADAM's docs/.vitepress/theme/index.ts, without its redirects): in the dark theme a page figure <img src=".../x.png"> shows x-dark.png, the silver-ink twin
// drawn by docs/scripts/figure_theme.py, when it exists (dark-figures.data.ts), and stays a washi card otherwise
// (custom.css). The verification catalog, cards and casts are styled in custom.css.
import DefaultTheme from 'vitepress/theme'
import { h, nextTick, watch } from 'vue'
import { onContentUpdated, useData } from 'vitepress'
import { data as darkFigures } from './dark-figures.data'
import './custom.css'

// The twins are swapped on the client only, after the content is mounted (the server-rendered page always carries the
// light figure, so hydration never mismatches): on every content update (first mount, route change) and theme change.
// The light source is kept in data-light-src; a twin that fails to load falls back to the washi card.
const twins = new Set<string>(darkFigures)
const TWIN = 'seth-dark-twin'

function swapFigures(dark: boolean, base: string) {
  if (typeof document === 'undefined') return
  document.querySelectorAll<HTMLImageElement>('.vp-doc img').forEach((img) => {
    const light = img.dataset.lightSrc ?? img.getAttribute('src') ?? ''
    if (!img.dataset.lightSrc) {
      const path = new URL(light, window.location.href).pathname
      if (!path.startsWith(base) || !twins.has(path.slice(base.length))) return
      img.dataset.lightSrc = light
    }
    if (dark && !img.classList.contains(TWIN)) {
      img.addEventListener('error', () => {
        img.classList.remove(TWIN)
        img.setAttribute('src', light)
      }, { once: true })
      img.classList.add(TWIN)
      img.setAttribute('src', light.replace(/\.png$/, '-dark.png'))
    } else if (!dark && img.classList.contains(TWIN)) {
      img.classList.remove(TWIN)
      img.setAttribute('src', light)
    }
  })
}

const Layout = {
  setup() {
    const { site, isDark } = useData()
    onContentUpdated(() => swapFigures(isDark.value, site.value.base))
    watch(isDark, (dark) => nextTick(() => swapFigures(dark, site.value.base)))
    return () => h(DefaultTheme.Layout)
  },
}

export default { extends: DefaultTheme, Layout }
