import * as mat4 from './toji-gl-matrix-1f872b8/src/mat4.js'
function multiplyMatrices(m1, m2) {
    const result = new Array(9);
    for (let i = 0; i < 3; i++) {
        for (let j = 0; j < 3; j++) {
            // result[i*3 + j] = m1[i*3]*m2[j] + m1[i*3+1]*m2[3+j] + m1[i*3+2]*m2[6+j]
            result[i * 3 + j] = m1[i * 3] * m2[j] + 
                                m1[i * 3 + 1] * m2[3 + j] + 
                                m1[i * 3 + 2] * m2[6 + j];
        }
    }
    return result;
}
// p: vec3 of the point to test (e.g., sphere center)
// box: { center: vec3, rotation: 3x3 matrix (array of arrays), size: vec3 }

export function transformToLocal(p, box) {
    // Subtract box center
    const dx = p[0] - box.pos[0];
    const dy = p[1] - box.pos[1];
    const dz = p[2] - box.pos[2];

    // Apply transpose of rotation matrix (inverse rotation)
    const R = box.rot;
    return [
        dx * R[0][0] + dy * R[1][0] + dz * R[2][0],
        dx * R[0][1] + dy * R[1][1] + dz * R[2][1],
        dx * R[0][2] + dy * R[1][2] + dz * R[2][2],
    ];
}
function clamp(val, min, max) {
    return Math.max(min, Math.min(max, val));
}

export function closestPointOnBox(localPos, box) {
    const hx = box.sz[0] / 2;
    const hy = box.sz[1] / 2;
    const hz = box.sz[2] / 2;

    return [
        clamp(localPos[0], -hx, hx),
        clamp(localPos[1], -hy, hy),
        clamp(localPos[2], -hz, hz)
    ];
}
function dot(a,b) {
    return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
}

function cross(a,b) {
    return [
        a[1]*b[2]-a[2]*b[1],
        a[2]*b[0]-a[0]*b[2],
        a[0]*b[1]-a[1]*b[0]
    ];
}

function projectBoxOntoAxis(box, axis) {
    let r = 0;
    const axes = box.rot;
    for (let i=0;i<3;i++) {
        r += box.halfExtents[i] * Math.abs(dot(axes[i], axis));
    }
    const c = dot(box.pos, axis);
    return [c - r, c + r];
}

function overlap(proj1, proj2) {
    return proj1[0] <= proj2[1] && proj2[0] <= proj1[1];
}
function overlapAmount(proj1, proj2) {
    // proj = [min,max]
    return Math.min(proj1[1], proj2[1]) - Math.max(proj1[0], proj2[0]);
}

function obbCollision(boxA, boxB) {
    const axesA = boxA.rot;
    const axesB = boxB.rot;

    const axes = [];
    // 3 axes from A
    axes.push(...axesA);
    // 3 axes from B
    axes.push(...axesB);
    // 9 cross products
    for(let i=0;i<3;i++){
        for(let j=0;j<3;j++){
            const cp = cross(axesA[i], axesB[j]);
            // skip near-zero vectors
            if (cp[0]!==0 || cp[1]!==0 || cp[2]!==0)
                axes.push(cp);
        }
    }

    // test all axes
    for (let axis of axes) {
        const projA = projectBoxOntoAxis(boxA, axis);
        const projB = projectBoxOntoAxis(boxB, axis);
        if (!overlap(projA, projB)) return false; // separating axis found
    }
    return true; // no separating axis → collision
}
function obbCollisionWithMTV(boxA, boxB) {
    const axesA = boxA.rot;
    const axesB = boxB.rot;

    const axes = [];
    axes.push(...axesA);
    axes.push(...axesB);

    for(let i=0;i<3;i++){
        for(let j=0;j<3;j++){
            const cp = cross(axesA[i], axesB[j]);
            if (cp[0]!==0 || cp[1]!==0 || cp[2]!==0)
                axes.push(cp);
        }
    }

    let minOverlap = Infinity;
    let mtvAxis = null;

    for (let axis of axes) {
        // normalize axis
        const len = Math.hypot(...axis);
        if(len === 0) continue;
        axis = axis.map(v => v/len);

        const projA = projectBoxOntoAxis(boxA, axis);
        const projB = projectBoxOntoAxis(boxB, axis);
        if (!overlap(projA, projB)) return null; // no collision

        const o = overlapAmount(projA, projB);
        if (o < minOverlap) {
            minOverlap = o;
            mtvAxis = axis;
        }
    }

    return { axis: mtvAxis, overlap: minOverlap };
}
export function sepBoxes(boxA,boxB){
    const mtv = obbCollisionWithMTV(boxA, boxB);
if (mtv) {
    const direction = [
        boxB.pos[0] - boxA.pos[0],
        boxB.pos[1] - boxA.pos[1],
        boxB.pos[2] - boxA.pos[2]
    ];

    // Ensure axis points from A to B
    const dotProd = direction[0]*mtv.axis[0] + direction[1]*mtv.axis[1] + direction[2]*mtv.axis[2];
    const sepAxis = dotProd < 0 ? mtv.axis.map(v => -v) : mtv.axis;

    // Move boxes apart
    // For example, move B fully, or half each
    const halfOverlap = mtv.overlap / 2;

    boxA.pos = boxA.pos.map((v,i) => v - sepAxis[i]*halfOverlap);
    boxB.pos = boxB.pos.map((v,i) => v + sepAxis[i]*halfOverlap);
}

}
export function eulerAngles(x, y, z) {
    const xm = [
        1, 0, 0,
        0, Math.cos(x), -Math.sin(x),
        0, Math.sin(x), Math.cos(x)
    ];
    const ym = [
        Math.cos(y), 0, Math.sin(y),
        0, 1, 0,
        -Math.sin(y), 0, Math.cos(y)
    ];
    const zm = [
        Math.cos(z), -Math.sin(z), 0,
        Math.sin(z), Math.cos(z), 0,
        0, 0, 1
    ];

    // Combine matrices correctly using matrix multiplication.
    // The order depends on the desired Euler angle sequence (e.g., Z * Y * X)
    // This example uses Z * Y * X order (common for extrinsic rotations):
    let rm = multiplyMatrices(zm, ym);
    rm = multiplyMatrices(rm, xm);
    
    return rm;
}

export function multMat(v,m){
    const x = m[0]*v.x+m[1]*v.y+m[2]*v.z
    const y = m[3]*v.x+m[4]*v.y+m[5]*v.z
    const z = m[6]*v.x+m[7]*v.y+m[8]*v.z
    return {x,y,z}
}
export class triangleBuffer {
    inds = []
    verts = []
    vbuffer = []
    ibuffer = []
    gl
    ind = 0
    boxes=[]
    balls=[]
    /**
     * 
     * @param {WebGL2RenderingContext} gl 
     */
    constructor(gl) {
        this.vbuffer = null
        this.ibuffer = null
        this.gl = gl
        this.indexCount = 0
    }
    updateBuffers(){
        // bind provided buffers if available
        if(this.vbo) this.gl.bindBuffer(this.gl.ARRAY_BUFFER, this.vbo)
        if(this.ibo) this.gl.bindBuffer(this.gl.ELEMENT_ARRAY_BUFFER, this.ibo)
        this.gl.bufferData(this.gl.ARRAY_BUFFER, new Float32Array(this.verts), this.gl.DYNAMIC_DRAW)
        this.gl.bufferData(this.gl.ELEMENT_ARRAY_BUFFER, new Uint16Array(this.inds), this.gl.DYNAMIC_DRAW)
        this.indexCount = this.inds.length
    }
    addQuad(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,nx,ny,nz,r,g,b){
        this.verts.push(
            x1,y1,z1,r,g,b,nx,ny,nz,
            x2,y2,z2,r,g,b,nx,ny,nz,
            x3,y3,z3,r,g,b,nx,ny,nz,
            x4,y4,z4,r,g,b,nx,ny,nz
        )
        this.inds.push(this.ind,this.ind+1,this.ind+2,this.ind,this.ind+2,this.ind+3)
        this.ind+=4
    }
    addBox(x,y,z,w,h,d,color){
        const x0 = x - w/2
        const x1 = x + w/2
        const y0 = y - h/2
        const y1 = y + h/2
        const z0 = z - d/2
        const z1 = z + d/2
        const baseind = this.ind // eh 1 character shorter
        const s3 = 1/Math.sqrt(3)
            const [r,g,b] = [color.r, color.g, color.b];

            this.addQuad(
                x0, y0, z1,
                x1, y0, z1,
                x1, y1, z1,
                x0, y1, z1,
                0, 0, 1,
                r, g, b
            );

            this.addQuad(
                x1, y0, z0,
                x0, y0, z0,
                x0, y1, z0,
                x1, y1, z0,
                0, 0, -1,
                r, g, b
            );

            this.addQuad(
                x0, y0, z0,
                x0, y0, z1,
                x0, y1, z1,
                x0, y1, z0,
                -1, 0, 0,
                r, g, b
            );

            this.addQuad(
                x1, y0, z1,
                x1, y0, z0,
                x1, y1, z0,
                x1, y1, z1,
                1, 0, 0,
                r, g, b
            );

            this.addQuad(
                x0, y1, z1,
                x1, y1, z1,
                x1, y1, z0,
                x0, y1, z0,
                0, 1, 0,
                r, g, b
            );

            this.addQuad(
                x0, y0, z0,
                x1, y0, z0,
                x1, y0, z1,
                x0, y0, z1,
                0, -1, 0,
                r, g, b
            );
            this.boxes.push({start:this.ind-24,end:this.ind,sverts:this.verts.slice((this.ind-24)*9,this.ind*9)})
        this.updateBuffers()
    }
    moveVerts(startInd,endInd,x,y,z){
        for (let i = startInd; i < endInd; i++){
            const vi = i*9
            this.verts[vi]+=x
            this.verts[vi+1]+=y
            this.verts[vi+2]+=z
        }
    }
    moveBoxTo(boxn,x,y,z){
        const box = this.boxes[boxn]
        for (let vi = box.start*9; vi < box.end*9; vi+=9){
            const i = vi-box.start*9
            const vert = box.sverts.slice(i,i+9)
            this.verts[vi] = vert[0]+x
            this.verts[vi+1] = vert[1]+y
            this.verts[vi+2] = vert[2]+z
        }
        // removed this.updateBuffers() so caller can batch uploads once per frame
    }
    rotBoxTo(boxn,matrix,x,y,z){
        /*moveBoxTo(boxn,x,y,z){
        const box = this.boxes[boxn]
        for (let vi = box.start*9; vi < box.end*9; vi+=9){
            const i = vi-box.start*9
            const vert = box.sverts.slice(i,i+9)
            this.verts[vi] = vert[0]+x
            this.verts[vi+1] = vert[1]+y
            this.verts[vi+2] = vert[2]+z
        }
        // removed this.updateBuffers() so caller can batch uploads once per frame
    } */
        const box = this.boxes[boxn]
        for (let vi = box.start*9; vi < box.end*9; vi+=9){
            const i = vi-box.start*9
            const vert = box.sverts.slice(i,i+9)
            const np = multMat({x:vert[0],y:vert[1],z:vert[2]},matrix)
            const nn = multMat({x:vert[6],y:vert[7],z:vert[8]},matrix)
            this.verts[vi] = np.x+x
            this.verts[vi+1] = np.y+y
            this.verts[vi+2] = np.z+z
            this.verts[vi+6] = nn.x
            this.verts[vi+7] = nn.y
            this.verts[vi+8] = nn.z
        }
    }
    rotBallTo(balln,matrix,x,y,z){
        const ball = this.balls[balln]
        for (let vi = ball.start*9; vi < ball.end*9; vi+=9){
            const i = vi-ball.start*9
            const vert = ball.sverts.slice(i,i+9)
            const np = multMat({x:vert[0],y:vert[1],z:vert[2]},matrix)
            const nn = multMat({x:vert[6],y:vert[7],z:vert[8]},matrix)
            this.verts[vi] = np.x+x
            this.verts[vi+1] = np.y+y
            this.verts[vi+2] = np.z+z
            this.verts[vi+6] = nn.x
            this.verts[vi+7] = nn.y
            this.verts[vi+8] = nn.z
        }
    }

    addBall(rad,bx,by,bz,r,g,b){
        const verts = []
        const inds = []
        const numVerts = 9.00001
        let ind = this.ind
        verts.push(bx, 1+by, bz, r, g, b,0, 1, 0)

        // generate vertices per stack / slice
        for (let i = 0; i < numVerts - 1; i++)
        {
            const phi = Math.PI * (i + 1) / numVerts
            for (let j = 0; j < numVerts; j++)
            {
                const theta = 2.0 * Math.PI * j / numVerts
                const nx = Math.sin(phi) * Math.cos(theta)
                const ny = Math.cos(phi)
                const nz = Math.sin(phi) * Math.sin(theta)
                const x = rad * nx+bx
                const y = rad * ny+by
                const z = rad * nz+bz
                verts.push(x, y, z, r, g, b, nx, ny, nz)
            }
        }

        // add bottom vertex
        verts.push(bx, by-1, bz, r, g, b, 0, -1, 0)

        // add top / bottom triangles
        
        for (let i = 0; i < numVerts - 1; i++)
        {
            let i0 = i + 1
            let i1 = (i + 1) % numVerts + 1
            inds.push(ind, i1+ind,i0+ind)
            i0 = i + numVerts * (numVerts - 2) + 1
            i1 = (i + 1) % numVerts + numVerts * (numVerts - 2) + 1
            inds.push(numVerts * (numVerts - 1) + 1+ind, i0+ind, i1+ind)
        }
        for (let j = 0; j < numVerts - 2; j++)
        {
            const j0 = j * numVerts + 1
            const j1 = (j + 1) * numVerts + 1
            for (let i = 0; i < numVerts; i++)
            {
                const i0 = j0 + i
                const i1 = j0 + (i + 1) % numVerts
                const i2 = j1 + (i + 1) % numVerts
                const i3 = j1 + i
                inds.push(i0+ind,i1+ind,i2+ind,i0+ind,i2+ind,i3+ind)
            }
        }
        this.verts.push(...verts)
        this.inds.push(...inds)
        const bl = verts.length/9
        this.ind = this.verts.length/9
        console.log(bl,ind,this.ind)
        this.balls.push({start:this.ind-bl,end:this.ind,sverts:this.verts.slice((this.ind-bl)*9,this.ind*9)})
        this.updateBuffers()
    }

    moveBallTo(balln,x,y,z){
        const ball = this.balls[balln]
        for (let vi = ball.start*9; vi < ball.end*9; vi+=9){
            const i = vi-ball.start*9
            const vert = ball.sverts.slice(i,i+9)
            this.verts[vi] = vert[0]+x
            this.verts[vi+1] = vert[1]+y
            this.verts[vi+2] = vert[2]+z
        }
        // removed this.updateBuffers() so caller can batch uploads once per frame
    }
}
export class Scene {
    buffer
    lightdir
    lightcol
    ambient
    camera
    gl
    lookat
    prog
    /**
     * 
     * @param {Array<Number>} lightdir 
     * @param {Array<Number>} lightcol 
     * @param {Array<Number>} ambient 
     * @param {Array<Number>} camera 
     * @param {Number} fov 
     * @param {WebGL2RenderingContext} gl 
     */
    constructor(lightdir,lightcol,ambient,camera,lookat,gl,prog){
        this.lightdir = lightdir
        this.lightcol = lightcol
        this.ambient = ambient
        this.camera=camera
        this.fov=fov
        this.gl=gl
        this.lookat=lookat
        this.prog=prog
        this.buffer = new triangleBuffer(gl)
        const ldloc = gl.getUniformLocation(prog, 'lightdir')
        const lcloc = gl.getUniformLocation(prog, 'lightColor')
        const acloc = gl.getUniformLocation(prog, 'ambientColor')
        gl.uniform3fv(ldloc,lightdir)
        gl.uniform3fv(lcloc,lightcol)
        gl.uniform3fv(acloc,ambient)
        const cameram = mat4.create()
        mat4.lookAt(camera, camera, lookat, [0, 1, 0])
        const vmatloc = gl.getUniformLocation(prog, 'vmat')
        gl.uniformMatrix4fv(vmatloc, false, cameram)
    }

}
export class Camera {
    pos = [0,5,-10]
    lookat = [0,0,0]
    gl
    prog
    cammat
    lookdir=[0,0,0]
    /**
     * 
     * @param {Array<Number>} pos 
     * @param {Array<Number>} lookingat 
     * @param {WebGL2RenderingContext} gl 
     */
    constructor(pos,lookingat,gl,program){
        this.pos = pos
        this.lookat = lookingat
        this.gl=gl
        this.prog=program
        this.cammat = mat4.create()
        this.getLookDir()
        mat4.lookAt(this.cammat,pos,lookingat,[0,1,0])
    }
    getLookDir(){
        const dist = Math.sqrt((this.pos[0]-this.lookat[0])**2+(this.pos[1]-this.lookat[1])**2+(this.pos[2]-this.lookat[2])**2)

        this.lookdir[0] = (this.pos[0]-this.lookat[0])/dist
        this.lookdir[1] = (this.pos[1]-this.lookat[1])/dist
        this.lookdir[2] = (this.pos[2]-this.lookat[2])/dist
    }
    updateCam(){

        mat4.lookAt(this.cammat,this.pos,this.lookat,[0,1,0])
        const vmatloc = this.gl.getUniformLocation(this.prog, 'vmat')
        this.gl.uniformMatrix4fv(vmatloc, false, this.cammat)
        this.getLookDir()
    }
    moveCamTo(pos){
        this.pos = pos
        this.getLookDir()
    }
    moveCam(movement){
        this.pos = [this.pos[0]+movement[0],this.pos[1]+movement[1],this.pos[2]+movement[2]]
        this.lookat = [this.lookat[0]+movement[0],this.lookat[1]+movement[1],this.lookat[2]+movement[2]]
        this.getLookDir()
    }
    moveCam2(movement){
        // Interpret movement as [forward, up, strafe]
        // rotate the forward/strafing to match the camera's yaw (horizontal look direction)
        
            const forwardInput = movement[2] || 0
            const upInput = movement[1] || 0
            const strafeInput = movement[0] || 0

            // world-space forward = lookat - pos, flattened to XZ to ignore pitch for horizontal motion
            let fx = this.lookat[0] - this.pos[0]
            let fz = this.lookat[2] - this.pos[2]
            const fh = Math.hypot(fx, fz)
            if (fh > 1e-6) {
                fx /= fh
                fz /= fh
            } else {
                fx = 0
                fz = 1
            }

            // right = normalize(cross(forward, up))
            // with up = [0,1,0], right = [fz, 0, -fx]
            let rx = fz
            let ry = 0
            let rz = -fx
            const rd = Math.hypot(rx, ry, rz)
            if (rd > 1e-6) {
                rx /= rd
                ry /= rd
                rz /= rd
            } else {
                rx = 1; ry = 0; rz = 0
            }

            // compose world movement
            movement = [
                fx * forwardInput + rx * strafeInput,
                upInput,
                fz * forwardInput + rz * strafeInput
            ]
        

        this.pos = [this.pos[0]+movement[0],this.pos[1]+movement[1],this.pos[2]+movement[2]]
        this.lookat = [this.lookat[0]+movement[0],this.lookat[1]+movement[1],this.lookat[2]+movement[2]]
        this.getLookDir()

    }
    rotCam(angle){
        try{
        this.lookdir = [Math.sin(angle),0,Math.cos(angle)]
        this.lookat = [this.pos[0]+this.lookdir[0],this.pos[1]+this.lookdir[1],this.pos[2]+this.lookdir[2]]
       
        }catch(e){alert(e)}
    }
}