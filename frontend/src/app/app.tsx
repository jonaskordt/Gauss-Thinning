import React, { useCallback, useEffect, useState } from "react";
import styled from "styled-components";
import { GlobalStyles } from "../theme";
import { Renderer } from "../rendering";
import { ImportZone } from "../components";
import { useIsDraggedOver } from "../utils";

const renderer = new Renderer();

const Cover = styled.div`
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

  const [isObjectLoaded, setIsObjetLoaded] = useState(false);
  const importFile = useCallback((url: string) => {
    renderer.importObject(url);
    setIsObjetLoaded(true);
  }, []);

  const [isDraggedOver, { onDrop, ...dragListeners }] = useIsDraggedOver();
  const onOutsideDrop = useCallback(
    (event: React.DragEvent) => {
      event.preventDefault();
      onDrop();
    },
    [onDrop]
  );

  return (
    <Cover {...dragListeners} onDrop={onOutsideDrop}>
      <GlobalStyles />
      <Cover ref={setRef} />
      <ImportZone
        isDraggedOver={isDraggedOver}
        onDrop={onDrop}
        importFile={importFile}
        isObjectLoaded={isObjectLoaded}
      />
    </Cover>
  );
};
