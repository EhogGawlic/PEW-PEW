precision mediump float;
attribute vec3 aPosition;
attribute vec3 acolor;
attribute vec3 anorm;
varying vec3 vColor;
varying vec3 vNorm;
varying vec3 vPos;
uniform mat4 pmat;
uniform mat4 vmat;

void main(void) {
    gl_Position = pmat * vmat * vec4(aPosition, 1.0);
    vColor = acolor;
    vNorm = anorm;
    vPos = aPosition;
}