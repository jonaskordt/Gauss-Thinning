export class Client {
  protected webSocket?: WebSocket;
  protected isConnected = false;

  public connect() {
    this.webSocket = new WebSocket("ws://localhost:5678/");
    this.webSocket.onopen = () => {
      this.isConnected = true;
    };
    this.webSocket.onmessage = (event: MessageEvent) => console.log(event.data);
  }

  public sendObject(dataURL: string) {
    if (!this.webSocket || !this.isConnected) return;
    this.webSocket?.send(JSON.stringify({ object: dataURL }));
  }

  public sendPath(path: number[]) {
    if (!this.webSocket || !this.isConnected) return;
    this.webSocket?.send(JSON.stringify({ path }));
  }
}
