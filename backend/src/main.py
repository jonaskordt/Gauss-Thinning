import asyncio
import websockets


async def socket(websocket):
    while True:
        await websocket.send("test")
        await asyncio.sleep(3)


start_server = websockets.serve(socket, "localhost", 5678)

asyncio.get_event_loop().run_until_complete(start_server)
asyncio.get_event_loop().run_forever()
