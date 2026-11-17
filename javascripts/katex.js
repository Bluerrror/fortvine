// pymdownx.arithmatex (generic mode) wraps math as \( ... \) and \[ ... \]; render exactly those
document$.subscribe(({ body }) => {
  renderMathInElement(body, {
    delimiters: [
      { left: "\\(", right: "\\)", display: false },
      { left: "\\[", right: "\\]", display: true }
    ],
    throwOnError: false
  })
})
