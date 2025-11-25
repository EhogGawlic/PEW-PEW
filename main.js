import {Game} from "./teedee/engine.js"
import { Camera } from "./teedee/utils.js"
async function main(){
    const g = new Game()
    let t = 0
    g.init("canv",()=>{
        g.addBox(0,-10,0,1024,1,1024,g.things,"Baseplate",{r:1,g:1,b:1})
            g.things.Baseplate.moveTo(512,0,0)

        g.camera = new Camera([-1,4,10],[1,0,0],g.gl,g.prog)
        g.camera.updateCam()
    },
    ()=>{
        try{
        g.camera = new Camera([Math.sin(t*0.1)*10,4,Math.cos(t*0.1)*10],[1,0,0],g.gl,g.prog)
        g.camera.updateCam()
        t++

        }catch(e){alert(e)}
    })
}
main().catch(e=>{alert(e)})