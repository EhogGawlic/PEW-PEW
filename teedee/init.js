import { readTextFile } from "./filereader.js"
/**
 * 
 * @param {String} canvasid 
 * @returns {canvas: HTMLCanvasElement, gl: WebGL2RenderingContext}
 */
export function initCanvas(canvasid){
    const canvas = document.getElementById(canvasid)
    canvas.width = canvas.clientWidth
    canvas.height = canvas.clientHeight
    /**
     * @type {WebGLRenderingContext}
     */
    const gl = canvas.getContext('webgl2')
    
    gl.viewport(0, 0, canvas.width, canvas.height)
    gl.enable(gl.DEPTH_TEST)
    gl.enable(gl.CULL_FACE)
    gl.cullFace(gl.BACK)
    gl.frontFace(gl.CCW)
    return {canvas, gl}
}
/**
 * 
 * @param {WebGL2RenderingContext} gl 
 * @param {String} shaderType 
 * @returns {vs: WebGL2Shader, fs: WebGL2Shader}
 */
export async function initShaders(gl, shaderType){
    let vssp = './teedee/shaders/vsbasic.glsl'
    let fssp = './teedee/shaders/fsbasic.glsl'
    switch(shaderType){
        case 'shaded':
            vssp = './teedee/shaders/vsshaded.glsl'
            fssp = './teedee/shaders/fsshaded.glsl'
            break
        case 'tex':
            vssp = './teedee/shaders/vstex.glsl'
            fssp = './teedee/shaders/fstex.glsl'
            break
    }
    const vstxt = await readTextFile(vssp)
    const fstxt = await readTextFile(fssp)
    const vs = gl.createShader(gl.VERTEX_SHADER)
    const fs = gl.createShader(gl.FRAGMENT_SHADER)
    gl.shaderSource(vs, vstxt)
    gl.shaderSource(fs, fstxt)
    gl.compileShader(vs)
    gl.compileShader(fs)
    if(!gl.getShaderParameter(vs, gl.COMPILE_STATUS)){
        console.error('Vertex shader failed to compile')
        console.error(gl.getShaderInfoLog(vs))
    }
    if(!gl.getShaderParameter(fs, gl.COMPILE_STATUS)){
        console.error('Fragment shader failed to compile')
        console.error(gl.getShaderInfoLog(fs))
    }
    return {vs, fs}
}
/**
 * 
 * @param {WebGL2RenderingContext} gl 
 * @param {WebGL2Shader} vs 
 * @param {WebGL2Shader} fs 
 * @returns {WebGL2Program}
 */
export function createProgram(gl, vs, fs){
    console.log(vs,fs)
    const prog = gl.createProgram()
    gl.attachShader(prog, vs)
    gl.attachShader(prog, fs)
    gl.linkProgram(prog)
    gl.validateProgram(prog)
    if(!gl.getProgramParameter(prog, gl.LINK_STATUS)){
        console.error('Program failed to link')
        console.error(gl.getProgramInfoLog(prog))
    }
    if(!gl.getProgramParameter(prog, gl.VALIDATE_STATUS)){
        console.error('Program failed to validate')
        console.error(gl.getProgramInfoLog(prog))
    }
    return prog
}
