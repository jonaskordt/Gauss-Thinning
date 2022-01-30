import { createGlobalStyle } from "styled-components";

export const GlobalStyles = createGlobalStyle`
  body {
    font-family: "Gothic A1", -apple-system, BlinkMacSystemFont, "Segoe UI",
      "Roboto", "Oxygen", "Ubuntu", "Cantarell", "Fira Sans", "Droid Sans",
      "Helvetica Neue", sans-serif;
    -webkit-font-smoothing: antialiased;
    -moz-osx-font-smoothing: grayscale;
    height: 100%;
    margin: 0;
  }

  html {
    height: 100%;
    min-height: 100vh;
  }

  main,
  p,
  h1,
  h2,
  h3,
  h4,
  h5,
  h6,
  hr {
    margin: 0;
  }

  input {
    border-radius: 0;
  }

  #root {
    background-color: #FFFFFF;
    color: #000000;
    height: 100%;
    overflow: auto;
  }

`;
