import asyncio
from typing import Optional
from urllib.request import urlopen
import websockets
import json
import igl
import os
import numpy as np

v: Optional[np.array] = None
f: Optional[np.array] = None


async def socket(websocket):
    async for message in websocket:
        data = json.loads(message)
        if "object" in data:
            content = urlopen(data["object"])
            tmp_file_name = "tmp.obj"
            file = open(tmp_file_name, "w")
            file.write(content.read().decode("utf-8"))
            file.close()
            v, f = igl.read_triangle_mesh(tmp_file_name)
            os.remove(tmp_file_name)
        if "path" in data:
            if v is None or f is None:
                return
            path = data["path"]
            print(path)

            # await websocket.send("test")


start_server = websockets.serve(socket, "localhost", 5678)

asyncio.get_event_loop().run_until_complete(start_server)
asyncio.get_event_loop().run_forever()
