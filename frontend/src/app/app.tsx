import React, { useEffect, useState } from "react";
import styled from "styled-components";
import { GlobalStyles } from "../theme";
import { Renderer } from "../rendering";

const renderer = new Renderer();

const CanvasContainer = styled.div`
  position: absolute;
  top: 0;
  bottom: 0;
  right: 0;
  left: 0;
`;

export const App: React.FC = () => {
  const [ref, setRef] = useState<HTMLDivElement | null>(null);
  const canvas = renderer.renderer.domElement;
  useEffect(() => {
    if (ref && canvas) {
      ref.appendChild(canvas);
      renderer.resize();
    }

    return () => {
      if (ref) ref.innerHTML = "";
    };
  }, [canvas, ref]);

  return (
    <>
      <GlobalStyles />
      <CanvasContainer ref={setRef} />
    </>
  );
};
