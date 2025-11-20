precision mediump float;
attribute vec3 aPosition;
attribute vec3 acolor;
varying vec3 vColor;

uniform mat4 pmat;
uniform mat4 vmat;

void main(void) {
    gl_Position = pmat * vmat * vec4(aPosition, 1.0);
    vColor = acolor;
}