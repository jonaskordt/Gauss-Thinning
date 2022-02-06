import * as THREE from "three";
import { OrbitControls } from "three/examples/jsm/controls/OrbitControls";
import { OBJLoader } from "three/examples/jsm/loaders/OBJLoader";

export class Renderer {
  public renderer = new THREE.WebGLRenderer({ antialias: true, alpha: true });

  protected scene = new THREE.Scene();
  protected camera = new THREE.PerspectiveCamera(
    60,
    window.innerWidth / window.innerHeight,
    0.0001,
    10
  );

  protected orbitControls: OrbitControls;
  protected object?: THREE.Object3D;

  protected loader = new OBJLoader();

  protected lazyRenderTriggered = true;

  constructor() {
    const directionalLight = new THREE.DirectionalLight("white", 0.6);
    directionalLight.position.y = 1;
    this.camera.add(directionalLight);
    const ambientLight = new THREE.AmbientLight("white", 0.5);
    this.scene.add(this.camera, ambientLight);

    this.camera.position.set(2, 2, 2);
    this.camera.lookAt(0, 0, 0);

    this.orbitControls = new OrbitControls(
      this.camera,
      this.renderer.domElement
    );
    this.orbitControls.minDistance = 0.5;
    this.orbitControls.maxDistance = 8;
    this.orbitControls.addEventListener("change", this.lazyRender);

    window.addEventListener("resize", this.resize);
    this.resize();

    this.renderer.setAnimationLoop(this.animate);
  }

  public resize = () => {
    const aspect = window.innerWidth / window.innerHeight;

    this.camera.aspect = aspect;
    this.camera.updateProjectionMatrix();

    this.renderer.setSize(window.innerWidth, window.innerHeight);

    this.lazyRender();
  };

  public lazyRender = () => {
    this.lazyRenderTriggered = true;
  };

  private eagerRender = () => {
    this.renderer.render(this.scene, this.camera);
  };

  private animate = () => {
    if (this.lazyRenderTriggered) {
      this.lazyRenderTriggered = false;
      this.eagerRender();
    }
  };

  public loadObject = (url: string) => {
    this.loader.load(url, (object) => {
      if (this.object) {
        this.scene.remove(this.object);
      }

      this.object = object;
      this.scene.add(object);

      this.lazyRender();
    });
  };
}
