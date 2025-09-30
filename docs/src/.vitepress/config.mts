import { defineConfig } from 'vitepress'
import { tabsMarkdownPlugin } from 'vitepress-plugin-tabs'
import mathjax3 from "markdown-it-mathjax3";
import footnote from "markdown-it-footnote";
import path from 'path'

function getBaseRepository(base: string): string {
  if (!base || base === '/') return '/';
  const parts = base.split('/').filter(Boolean);
  return parts.length > 0 ? `/${parts[0]}/` : '/';
}

const baseTemp = {
  base: '/',// TODO: replace this in makedocs!
}

const navTemp = {
  nav: [
{ text: 'Home', link: '/index' },
{ text: 'Getting Started', collapsed: false, items: [
{ text: 'Julia 101', link: '/intro/intro_julia' },
{ text: 'Porosity.jl : Introduction', link: '/intro/getting_started' },
{ text: 'Visualization', link: '/intro/intro_figs' },
{ text: 'Structure', link: '/intro/structure' }]
 },
{ text: 'Models', collapsed: false, items: [
{ text: 'Conductivity models', link: '/models/conductivity' },
{ text: 'Elasticity models', link: '/models/elasticity' },
{ text: 'Viscosity models', link: '/models/viscosity' },
{ text: 'Anelasticity models', link: '/models/anelasticity' },
{ text: 'Adding a new method', link: '/models/new_method' }]
 },
{ text: 'Tutorials', collapsed: false, items: [
{ text: 'Mixing phases', link: '/tutorials/mixing_phases' },
{ text: 'Multi rock physics', link: '/tutorials/combine_models' },
{ text: 'Tuning rock physics hyperparameters', link: '/tutorials/tune_rp' },
{ text: 'Stochastic inversion', link: '/tutorials/stochastic_inverse' },
{ text: 'Automatic Differentiation', link: '/tutorials/ad' }]
 },
{ text: 'Solidus', link: '/solidus' },
{ text: 'API', link: '/api' }
]
,
}

const nav = [
  ...navTemp.nav,
  {
    component: 'VersionPicker'
  }
]

// https://vitepress.dev/reference/site-config
export default defineConfig({
  base: '/',// TODO: replace this in makedocs!
  title: 'Porosity.jl',
  description: 'Documentation for Porosity.jl',
  lastUpdated: true,
  cleanUrls: true,
  outDir: '../1', // This is required for MarkdownVitepress to work correctly...
  head: [
    
    ['script', {src: `${getBaseRepository(baseTemp.base)}versions.js`}],
    // ['script', {src: '/versions.js'], for custom domains, I guess if deploy_url is available.
    ['script', {src: `${baseTemp.base}siteinfo.js`}]
  ],
  
  vite: {
    define: {
      __DEPLOY_ABSPATH__: JSON.stringify('/'),
    },
    resolve: {
      alias: {
        '@': path.resolve(__dirname, '../components')
      }
    },
    optimizeDeps: {
      exclude: [ 
        '@nolebase/vitepress-plugin-enhanced-readabilities/client',
        'vitepress',
        '@nolebase/ui',
      ], 
    }, 
    ssr: { 
      noExternal: [ 
        // If there are other packages that need to be processed by Vite, you can add them here.
        '@nolebase/vitepress-plugin-enhanced-readabilities',
        '@nolebase/ui',
      ], 
    },
  },
  markdown: {
    math: true,
    config(md) {
      md.use(tabsMarkdownPlugin),
      md.use(mathjax3),
      md.use(footnote)
    },
    theme: {
      light: "github-light",
      dark: "github-dark"}
  },
  themeConfig: {
    outline: 'deep',
    logo: { src: '/logo.png', width: 24, height: 24},
    search: {
      provider: 'local',
      options: {
        detailedView: true
      }
    },
    nav,
    sidebar: [
{ text: 'Home', link: '/index' },
{ text: 'Getting Started', collapsed: false, items: [
{ text: 'Julia 101', link: '/intro/intro_julia' },
{ text: 'Porosity.jl : Introduction', link: '/intro/getting_started' },
{ text: 'Visualization', link: '/intro/intro_figs' },
{ text: 'Structure', link: '/intro/structure' }]
 },
{ text: 'Models', collapsed: false, items: [
{ text: 'Conductivity models', link: '/models/conductivity' },
{ text: 'Elasticity models', link: '/models/elasticity' },
{ text: 'Viscosity models', link: '/models/viscosity' },
{ text: 'Anelasticity models', link: '/models/anelasticity' },
{ text: 'Adding a new method', link: '/models/new_method' }]
 },
{ text: 'Tutorials', collapsed: false, items: [
{ text: 'Mixing phases', link: '/tutorials/mixing_phases' },
{ text: 'Multi rock physics', link: '/tutorials/combine_models' },
{ text: 'Tuning rock physics hyperparameters', link: '/tutorials/tune_rp' },
{ text: 'Stochastic inversion', link: '/tutorials/stochastic_inverse' },
{ text: 'Automatic Differentiation', link: '/tutorials/ad' }]
 },
{ text: 'Solidus', link: '/solidus' },
{ text: 'API', link: '/api' }
]
,
    editLink: { pattern: "https://https://github.com/ayushinav/Porosity.jl/edit/main/docs/src/:path" },
    socialLinks: [
      { icon: 'github', link: 'https://github.com/ayushinav/Porosity.jl' }
    ],
    footer: {
      message: 'Made with <a href="https://luxdl.github.io/DocumenterVitepress.jl/dev/" target="_blank"><strong>DocumenterVitepress.jl</strong></a><br>',
      copyright: `© Copyright ${new Date().getUTCFullYear()}.`
    }
  }
})
