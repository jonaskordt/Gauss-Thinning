import React, { useCallback, useRef, useState } from "react";

export const useIsDraggedOver = () => {
  const [isDraggedOver, setIsDraggedOver] = useState(false);
  const dragTimerRef = useRef<ReturnType<typeof setTimeout>>();

  const onDragOver = useCallback(
    (event: React.DragEvent) => {
      event.preventDefault();
      setIsDraggedOver(true);
      if (dragTimerRef.current !== undefined) {
        clearTimeout(dragTimerRef.current);
        dragTimerRef.current = undefined;
      }
    },
    [setIsDraggedOver]
  );

  const onDragEnd = useCallback(() => {
    if (dragTimerRef.current !== undefined) {
      clearTimeout(dragTimerRef.current);
    }
    dragTimerRef.current = setTimeout(() => {
      setIsDraggedOver(false);
    }, 50);
  }, [setIsDraggedOver]);

  const onDrop = useCallback(() => {
    setIsDraggedOver(false);
    if (dragTimerRef.current !== undefined) {
      clearTimeout(dragTimerRef.current);
      dragTimerRef.current = undefined;
    }
  }, [setIsDraggedOver]);

  return [
    isDraggedOver,
    { onDragOver, onDragEnd, onDragLeave: onDragEnd, onDrop },
  ] as [
    boolean,
    {
      onDragOver: typeof onDragOver;
      onDragEnd: typeof onDragEnd;
      onDragLeave: typeof onDragEnd;
      onDrop: typeof onDragEnd;
    }
  ];
};
