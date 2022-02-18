import React, { useCallback, useEffect, useState } from "react";
import styled from "styled-components";
import { GlobalStyles } from "../theme";
import { ImportZone } from "../components";
import { useIsDraggedOver } from "../utils";
import { Client } from "../client";

const client = new Client();
const { renderer } = client;

const Cover = styled.div`
  position: absolute;
  top: 0;
  bottom: 0;
  right: 0;
  left: 0;
`;

const GlobalButton = styled.button`
  position: absolute;
  bottom: 20px;
  right: 20px;

  background: transparent;
  border: none;
  font-family: inherit;
  overflow: visible;
  text-transform: none;

  padding: 10px;
  background-color: lightgray;
  border-radius: 10px;
  font-size: 14px;
  border: 1px solid gray;
  cursor: pointer;
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
    client.importObject(url);
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

      {isObjectLoaded && (
        <GlobalButton onPointerDown={client.requestGlobalThinning}>
          Global Thinning
        </GlobalButton>
      )}
    </Cover>
  );
};
