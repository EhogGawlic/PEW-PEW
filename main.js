import {Game} from "./teedee/engine.js"
async function main(){
    const g = new Game()
    g.init("canv",()=>{
        g.addBox(0,0,0,1,1,1,g.things,"BOX",{r:1,g:0,b:1})
        
    })
}
main().catch(e=>{alert(e)})