import * as THREE from "three";
import { OrbitControls } from "three/examples/jsm/controls/OrbitControls";
import { OBJLoader } from "three/examples/jsm/loaders/OBJLoader";
import { Client } from "../client";

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
  protected object?: THREE.Mesh;
  protected material = new THREE.MeshStandardMaterial({
    vertexColors: true,
    side: THREE.DoubleSide,
  });

  protected raycaster = new THREE.Raycaster();
  protected path: number[] = [];
  protected isPointerDown = false;

  protected loader = new OBJLoader();

  protected lazyRenderTriggered = true;

  constructor(protected client: Client) {
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
    this.orbitControls.mouseButtons = {
      LEFT: -1,
      MIDDLE: THREE.MOUSE.ROTATE,
      RIGHT: THREE.MOUSE.ROTATE,
    };
    this.orbitControls.minDistance = 0.5;
    this.orbitControls.maxDistance = 8;
    this.orbitControls.addEventListener("change", this.lazyRender);

    window.addEventListener("resize", this.resize);
    this.resize();

    window.addEventListener("pointerdown", this.onPointerDown);
    window.addEventListener("pointermove", this.onPointerMove);
    window.addEventListener("pointerup", this.onPointerUp);

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
      if (!object.children.length) return;

      if (this.object) {
        this.scene.remove(this.object);
      }

      this.object = object.children[0] as THREE.Mesh;
      this.object.material = this.material;

      const count = this.object.geometry.attributes.position.count;
      this.object.geometry.setAttribute(
        "color",
        new THREE.BufferAttribute(new Float32Array(count * 3).fill(0.9), 3)
      );

      this.scene.add(this.object);

      this.lazyRender();
    });
  };

  public updateObject = (url: string) => {
    this.loader.load(url, (object) => {
      if (!this.object || !object.children.length) return;

      const update = object.children[0] as THREE.Mesh;
      const updatePositions = update.geometry.attributes.position;
      const updateNormals = update.geometry.attributes.normal;

      this.object.geometry.setAttribute("position", updatePositions);
      this.object.geometry.attributes.position.needsUpdate = true;
      this.object.geometry.setAttribute("normal", updateNormals);
      this.object.geometry.attributes.normal.needsUpdate = true;

      this.lazyRender();
    });
  };

  protected onPointerDown = (event: PointerEvent) => {
    if (event.button !== 0) return;

    this.isPointerDown = true;
    this.path = [];

    this.onPointerMove(event);
  };

  protected onPointerMove = (event: PointerEvent) => {
    if (!this.object || !this.isPointerDown) return;

    const x = (event.clientX / window.innerWidth) * 2 - 1;
    const y = -(event.clientY / window.innerHeight) * 2 + 1;
    this.raycaster.setFromCamera({ x, y }, this.camera);
    const intersections = this.raycaster.intersectObject(this.object);

    if (!intersections.length) return;

    const { face, faceIndex } = intersections[0];

    if (faceIndex === undefined || this.path.includes(faceIndex) || !face) {
      return;
    }

    [face.a, face.b, face.c].forEach((vertexIndex) => {
      this.object!.geometry.attributes.color.setXYZ(vertexIndex, 1, 0, 0);
    });
    this.object.geometry.attributes.color.needsUpdate = true;

    this.lazyRender();

    this.path.push(faceIndex);
  };

  protected onPointerUp = (event: PointerEvent) => {
    if (event.button !== 0) return;

    this.isPointerDown = false;

    if (!this.path.length) return;

    this.client.requestLocalThinning(this.path);
  };
}
