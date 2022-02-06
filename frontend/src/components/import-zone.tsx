import React, { useCallback } from "react";
import styled from "styled-components";

const Container = styled.div`
  position: absolute;
  top: 100px;
  bottom: 100px;
  right: 100px;
  left: 100px;
  border: 10px dashed #ccc;
  pointer-events: auto;
  display: flex;
  align-items: center;
  justify-content: center;
`;

const Text = styled.span`
  font-size: 30px;
  font-weight: 900;
  color: #ccc;
`;

const ImportNoteContainer = styled.div`
  position: absolute;
  top: 0;
  bottom: 0;
  right: 0;
  left: 0;
  display: flex;
  align-items: center;
  justify-content: center;
`;

export const ImportZone: React.FC<{
  isDraggedOver: boolean;
  onDrop: () => void;
  importFile: (url: string) => void;
  isObjectLoaded?: boolean;
}> = ({ isDraggedOver, onDrop, importFile, isObjectLoaded = false }) => {
  const drop = useCallback(
    (event: React.DragEvent<HTMLDivElement>) => {
      event.preventDefault();
      onDrop();

      const file = event.dataTransfer.items[0].getAsFile();
      if (!file || !file.name.endsWith(".obj")) return;
      importFile(URL.createObjectURL(file));
    },
    [onDrop, importFile]
  );

  return isDraggedOver ? (
    <Container onDrop={drop}>
      <Text>Drop .obj file here.</Text>
    </Container>
  ) : !isObjectLoaded ? (
    <ImportNoteContainer>
      <Text>Drag in an .obj file to import it.</Text>
    </ImportNoteContainer>
  ) : null;
};
