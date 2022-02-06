export class Client {
  public connect() {
    const ws = new WebSocket("ws://localhost:5678/");
    ws.onmessage = (event: MessageEvent) => console.log(event.data);
  }
}
