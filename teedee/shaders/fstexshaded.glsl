precision mediump float;
uniform vec3 lightColor;
uniform vec3 ambientColor;
varying vec3 vColor;
varying vec2 vTexcoord;
varying sampler2D vTexture;
varying vec3 vNorm;
void main(void) {
    
    vec4 texColor = texture2D(vTexture, vTexcoord);
    vColor = vColor * texColor.rgb;
    float diffuse = max(dot(normalize(vNorm), normalize(lightdir)), 0.0);

    gl_FragColor = vec4(vColor * diffuse lightColor + ambientColor, 1.0);
}