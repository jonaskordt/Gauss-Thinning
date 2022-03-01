import asyncio
from typing import List, Optional
from urllib.request import urlopen
import websockets
import json
import igl
import os
import numpy as np
import base64
from local_gauss_thinning import local_gauss_thinning, constraint_gauss_thinning
from gauss_thinning import center


# add custom pybind lib
import sys
sys.path.insert(0,'../third_party/DevelopableApproximationViaGaussImageThinning')
from mainParallel import gaussThinning

async def socket(websocket):
    v: Optional[np.array] = None
    f: Optional[np.array] = None
    touched_face_indices: Optional[List[int]] = None
    tmp_file_name = "tmp.obj"

    async def send_vertecies(v_dash, kind="update"):
        igl.write_obj(tmp_file_name, v_dash, f)
        prefix = "data:application/octet-stream;base64,"
        data = open(tmp_file_name, "rb")
        contents = data.read()
        data_url = prefix + base64.b64encode(contents).decode("utf-8")
        data.close()
        os.remove(tmp_file_name)
        message = json.dumps({"kind": kind, "object": data_url})
        await websocket.send(message)

    async for message in websocket:
        data = json.loads(message)
        if "request" not in data:
            return

        if data["request"] == "import":
            content = urlopen(data["object"])
            file = open(tmp_file_name, "w")
            file.write(content.read().decode("utf-8"))
            file.close()
            
            v, f = igl.read_triangle_mesh(tmp_file_name)
            touched_face_indices = []
            center(v)

            os.remove(tmp_file_name)
            await send_vertecies(v, kind="import")
        elif data["request"] == "local_thinning":
            if v is None or f is None:
                return
            path = data["path"]

            v, touched = await local_gauss_thinning(
                v, f, path=path, num_iterations=10, callback=send_vertecies
            )
            for face_id in touched:
                if face_id not in touched_face_indices:
                    touched_face_indices.append(face_id)
        elif data["request"] == "global_thinning":
            if v is None or f is None:
                return

            active_faces = (
                [
                    face_id
                    for face_id in range(len(f))
                    if face_id not in touched_face_indices
                ]
                if touched_face_indices is not None
                else range(len(f))
            )

            v2 = gaussThinning(v, f)
            await send_vertecies(v2)

start_server = websockets.serve(socket, "localhost", 5678, max_size=None)

asyncio.get_event_loop().run_until_complete(start_server)
asyncio.get_event_loop().run_forever()
