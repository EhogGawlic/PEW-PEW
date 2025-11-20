import * as graphics from "./teedee.js"
import * as utils from "./utils.js"
import * as physics from "./physics.js"
import { triangleBuffer, Scene, Camera } from "./utils.js"
class Shape{
    si = 0
    ei = 0
    pos=[0,0,0]
    sz=[1,1,1]
    ppos = [0,0,0]
    color={r:1,g:1,b:1}
    buf
    shape
    bi
    v = [0,0,0]
    anchored=true
    things
    /**
     * 
     * @param {String} shape 
     * @param {triangleBuffer} buffer 
     * @param {Array<Number>|Number} size 
     */
    constructor(shape,buffer,size,color,things){
        
        this.things=things
        this.sz=size
        this.si = buffer.ind+0
        this.buf=buffer
        this.shape=shape
        this.color=color
        switch(shape){
            case 'box':
                buffer.addBox(this.pos[0],this.pos[1],this.pos[2],this.sz[0],this.sz[1],this.sz[2],this.color)
                this.ei = buffer.ind+0
                this.bi = buffer.boxes.length-1
                break
            case 'ball':
                buffer.addBall(this.sz,this.pos[0],this.pos[1],this.pos[2],this.color.r,this.color.g,this.color.b)
                this.ei = buffer.ind+0
                this.bi =buffer.balls.length-1
        }
    }
    moveTo(x,y,z){
        switch(this.shape){
            case 'box':
                this.buf.moveBoxTo(this.bi,x,y,z)
                break
            case 'ball':
                this.buf.moveBallTo(this.bi,x,y,z)
        }
        this.buf.updateBuffers()
        this.pos = [x,y,z]
    }
    update(){
        if (!this.anchored){
            this.v = [
                this.pos[0]-this.ppos[0],
                this.pos[1]-this.ppos[1]-0.01,
                this.pos[2]-this.ppos[2],
            ]
            this.ppos = [...this.pos]
            this.pos =[
                this.pos[0]+this.v[0],
                this.pos[1]+this.v[1],
                this.pos[2]+this.v[2]
            ]
            this.moveTo(this.pos[0],this.pos[1],this.pos[2])
            //TEST COLLISIONS
            if (!this.anchored){
            if (this.shape == 'ball'){
                getDescendants(this.things).forEach(other=>{
                    if (other===this) return;
                    if (other.anchored){const EPS = 1e-4; // ignore tiny numerical penetrations
                                
                        if (other.shape=='box'){
                            const {closestPoint,minDist} = physics.closestPointOnBox(
                                this.pos,
                                this.sz,
                                other.pos,
                                other.sz
                            )
                            if (minDist<0){
                                // collision handling with tolerance and proper penetration correction
                                if (minDist < -EPS) {
                                    const moveVec = [
                                        this.pos[0] - closestPoint[0],
                                        this.pos[1] - closestPoint[1],
                                        this.pos[2] - closestPoint[2],
                                    ];
                                    const len = Math.hypot(moveVec[0], moveVec[1], moveVec[2]);
                                    // penetration depth (closestPoint routine returns negative when penetrating)
                                    const penDepth = -minDist;
                                    // collision normal fallback: up if vector is degenerate
                                    const normal = len > 1e-6 ? [moveVec[0]/len, moveVec[1]/len, moveVec[2]/len] : [0,1,0];

                                    // push out by penetration depth along normal
                                    this.pos[0] += normal[0] * penDepth;
                                    this.pos[1] += normal[1] * penDepth;
                                    this.pos[2] += normal[2] * penDepth;

                                    // reflect velocity along normal with restitution (0.0001)
                                    const restitution = 0.0001;
                                    const vn = this.v[0]*normal[0] + this.v[1]*normal[1] + this.v[2]*normal[2];
                                    // remove the normal component and add reversed damped component
                                    this.v[0] = this.v[0] - (1 + restitution) * vn * normal[0];
                                    this.v[1] = this.v[1] - (1 + restitution) * vn * normal[1];
                                    this.v[2] = this.v[2] - (1 + restitution) * vn * normal[2];
                                    this.ppos = [
                                        this.pos[0]-this.v[0],
                                        this.pos[1]-this.v[1],
                                        this.pos[2]-this.v[2]
                                    ]
                                    this.moveTo(this.pos[0],this.pos[1],this.pos[2])
                                }
                             }
                        }
                        if (other.shape=='ball'){

                            const dist = physics.vecDist(this.pos,other.pos)
                            const minDist = this.sz+other.sz
                            if (dist<minDist){
                                // simple elastic collision response
                                const normal = [
                                    (this.pos[0]-other.pos[0])/dist,
                                    (this.pos[1]-other.pos[1])/dist,
                                    (this.pos[2]-other.pos[2])/dist,
                                ]
                                normal[0] *= (dist-minDist)
                                normal[1] *= (dist-minDist)
                                normal[2] *= (dist-minDist)
                                // push balls apart
                                this.pos[0] -= normal[0]
                                this.pos[1] -= normal[1]
                                this.pos[2] -= normal[2]
                                this.moveTo(this.pos[0],this.pos[1],this.pos[2])
                            }
                        }
                    }

                    if (!other.anchored){const EPS = 1e-4; // ignore tiny numerical penetrations
                         
                        if (other.shape=='ball'){

                            const dist = physics.vecDist(this.pos,other.pos)
                            const minDist = this.sz+other.sz
                            if (dist<minDist){
                                // simple elastic collision response
                                const normal = [
                                    (this.pos[0]-other.pos[0])/dist,
                                    (this.pos[1]-other.pos[1])/dist,
                                    (this.pos[2]-other.pos[2])/dist,
                                ]
                                normal[0] *= (dist-minDist)
                                normal[1] *= (dist-minDist)
                                normal[2] *= (dist-minDist)
                                // push balls apart
                                this.pos[0] -= normal[0]/2
                                this.pos[1] -= normal[1]/2
                                this.pos[2] -= normal[2]/2
                                other.pos[0] += normal[0]/2
                                other.pos[1] += normal[1]/2
                                other.pos[2] += normal[2]/2
                                other.moveTo(other.pos[0],other.pos[1],other.pos[2])
                                this.moveTo(this.pos[0],this.pos[1],this.pos[2])
                            }
                        }
                    }
                })
            }
        }
            //collision?*/
        }
        
    }
}
function getDescendants(obj) {
  const result = [];
  const visited = new WeakSet();

  function traverse(currentObj) {
    if (!currentObj || typeof currentObj !== 'object') return;
    if (visited.has(currentObj)) return;
    visited.add(currentObj);

    for (const key in currentObj) {
      if (!Object.prototype.hasOwnProperty.call(currentObj, key)) continue;
      const value = currentObj[key];

      if (value instanceof Shape) {
        result.push(value);
        // Do NOT recurse into Shape instances (they hold a reference back to things)
        continue;
      }

      if (typeof value === 'object' && value !== null && !Array.isArray(value)) {
        traverse(value);
      }
    }
  }

  traverse(obj);
  return result;
}
export class Game{
    things={}
    scripts={}
    buffer
    camera
    gl
    constructor(){
    }
    async init(canvasid,onload){
        const {gl,vbo,ibo,stype,prog,buffer} = await graphics.initall(canvasid,"shaded")
        this.gl=gl
        this.vbo=vbo
        this.ibo=ibo
        this.prog = prog
        this.buffer=buffer
        this.camera = new Camera([0,30,1],[0,0,0],gl,prog)
        this.camera.updateCam()
        onload()
        this.render()
    }
    addBox(x,y,z,sx,sy,sz,parent,name,color){
        const box = new Shape("box",this.buffer,[sx,sy,sz],color,this.things)
        box.moveTo(x,y,z)
        box.ppos = [x,y,z]
        parent[name]=box
        return parent[name]
    }
    addBall(x,y,z,rad,parent,name,color){
        const ball = new Shape("ball",this.buffer,rad,color,this.things)
        ball.moveTo(x,y,z)
        ball.ppos = [x,y,z]
        parent[name]=ball
        return parent[name]
    }
    render = async()=>{
        const t0 = performance.now()
        getDescendants(this.things).forEach(v=>v.update())
        this.buffer.updateBuffers()
        await graphics.render(this.gl,this.prog,this.vbo,this.ibo,this.buffer)
        let rawDt = performance.now() - t0

        // compute raw FPS and smooth the FPS (0 < alpha < 1)
        const alpha = 0.08
        const rawFps = 1000 / rawDt
        if (this._smoothedFps === undefined) this._smoothedFps = rawFps
        this._smoothedFps = this._smoothedFps * (1 - alpha) + rawFps * alpha
        const fpsVal = this._smoothedFps

        fps.innerText = fpsVal.toFixed(1)
        requestAnimationFrame(this.render)
    }
    moveCam(x,y,z){

    }
}