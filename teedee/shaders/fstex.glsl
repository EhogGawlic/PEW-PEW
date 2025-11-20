precision mediump float;
uniform vec3 lightColor;
uniform vec3 ambientColor;
varying vec3 vColor;
varying vec2 vTexcoord;
uniform sampler2D uTexture;

void main(void) {
    vec4 texColor = texture2D(uTexture, vTexcoord);
    vec3 col = vColor * texColor.rgb;
    gl_FragColor = vec4(col * lightColor + ambientColor, texColor.a);
}