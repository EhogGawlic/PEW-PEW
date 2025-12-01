import {Game} from "./teedee/engine.js"
import { Camera, eulerAngles } from "./teedee/utils.js"
async function main(){
    const g = new Game()
    let t = 0
    g.init("canv",()=>{
        const b = g.addBox(0,-1,0,10,1,10,g.things,"Baseplate",{r:1,g:1,b:1})
        
        const bol = g.addBall(0,2,0,1,g.things,"BOL",{r:0,g:1,b:0})
        bol.anchored = false
        //ball.anchored = false
        g.camera = new Camera([0,4,10],[0,1,0],g.gl,g.prog)
        g.camera.updateCam()

    },
    ()=>{
        try{
        //g.camera = new Camera([Math.sin(t*0.1)*10,4,Math.cos(t*0.1)*10],[1,0,0],g.gl,g.prog)
        //g.camera.updateCam()
        t++
        }catch(e){alert(e)}
    })
}
main().catch(e=>{alert(e)})