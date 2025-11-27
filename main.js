import {Game} from "./teedee/engine.js"
import { Camera, eulerAngles } from "./teedee/utils.js"
async function main(){
    const g = new Game()
    let t = 0
    g.init("canv",()=>{
       // const b = g.addBox(0,0,0,10,1,10,g.things,"Baseplate",{r:1,g:1,b:1})
        //b.rotTo(eulerAngles(1,0,0))
        const ball = g.addBall(0,50,0,1,g.things,"BOL",{r:1,g:0,b:0})
        ball.anchored = false
        g.camera = new Camera([0,4,10],[0,1,0],g.gl,g.prog)
        g.camera.updateCam()
    },
    ()=>{
        try{
        g.camera = new Camera([Math.sin(t*0.1)*10,4,Math.cos(t*0.1)*10],[1,0,0],g.gl,g.prog)
        //g.camera.updateCam()
        //g.things.Baseplate.rotTo(eulerAngles(0,t*5.9632,0))
        t++
            console.log(g.things.BOL.pos)
        }catch(e){alert(e)}
    })
}
main().catch(e=>{alert(e)})