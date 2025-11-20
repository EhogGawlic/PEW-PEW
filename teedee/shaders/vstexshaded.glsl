precision mediump float;
attribute vec3 aPosition;
attribute vec3 acolor;
attribute vec2 atexcoord;
attribute vec3 anorm;
attribute sampler2D aTexture;
varying vec3 vColor;
varying vec2 vTexcoord;
varying sampler2D vTexture;
varying vec3 vNorm;
uniform mat4 pmat;
uniform mat4 vmat;

void main(void) {
    gl_Position = pmat * vmat * vec4(aPosition, 1.0);
    vColor = acolor;
    vTexcoord = atexcoord;
    vTexture = aTexture;
    vNorm = anorm;
}