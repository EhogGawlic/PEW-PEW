precision mediump float;
uniform vec3 lightpos;
uniform vec3 lightColor;
uniform vec3 ambientColor;
varying vec3 vColor;
varying vec3 vNorm;
varying vec3 vPos;
void main(void) {
    float diffuse = max(dot(normalize(vNorm), normalize(vec3(lightpos)-vPos)), 0.0);
    gl_FragColor = vec4(vColor * diffuse*lightColor + ambientColor, 1.0);
}