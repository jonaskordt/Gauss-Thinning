import { Renderer } from "../rendering";

export class Client {
  protected webSocket?: WebSocket;
  protected isConnected = false;

  public renderer: Renderer;

  constructor() {
    this.connect();
    this.renderer = new Renderer(this);
  }

  protected connect = () => {
    this.webSocket = new WebSocket("ws://localhost:5678/");
    this.webSocket.onopen = () => {
      this.isConnected = true;
    };
    this.webSocket.onerror = () => setTimeout(this.connect, 1000);
    this.webSocket.onmessage = this.handleMessage;
  };

  public importObject(dataURL: string) {
    if (!this.webSocket || !this.isConnected) return;
    this.webSocket?.send(
      JSON.stringify({ request: "import", object: dataURL })
    );
  }

  public requestLocalThinning = (path: number[]) => {
    if (!this.webSocket || !this.isConnected) return;
    this.webSocket?.send(JSON.stringify({ request: "local_thinning", path }));
  };

  public requestGlobalThinning = () => {
    if (!this.webSocket || !this.isConnected) return;
    this.webSocket?.send(JSON.stringify({ request: "global_thinning" }));
  };

  private handleMessage = (event: MessageEvent) => {
    const message = JSON.parse(event.data);

    const data_url = message.object;

    switch (message.kind) {
      case "import":
        this.renderer.loadObject(data_url);
        return;
      case "update":
        this.renderer.updateObject(data_url);
        return;
      default:
        return;
    }
  };
}
