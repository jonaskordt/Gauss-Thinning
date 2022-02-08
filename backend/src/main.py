import asyncio
from typing import Optional
from urllib.request import urlopen
import websockets
import json
import igl
import os
import numpy as np
import base64
from local_gauss_thinning import local_gauss_thinning

v: Optional[np.array] = None
f: Optional[np.array] = None
tmp_file_name = "tmp.obj"


async def socket(websocket):
    async def send_updated_vertecies(v_dash):
        igl.write_obj(tmp_file_name, v_dash, f)
        prefix = "data:application/octet-stream;base64,"
        data = open(tmp_file_name, "rb")
        contents = data.read()
        data_url = prefix + base64.b64encode(contents).decode("utf-8")
        data.close()
        os.remove(tmp_file_name)
        message = json.dumps({"update": data_url})
        await websocket.send(message)

    async for message in websocket:
        data = json.loads(message)
        if "object" in data:
            content = urlopen(data["object"])
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

            v = await local_gauss_thinning(
                v, f, path=path, num_iterations=15, callback=send_updated_vertecies
            )


start_server = websockets.serve(socket, "localhost", 5678)

asyncio.get_event_loop().run_until_complete(start_server)
asyncio.get_event_loop().run_forever()
