export function vecDist(a,b){
    return Math.sqrt(
        (a[0]-b[0])**2+
        (a[1]-b[1])**2+
        (a[2]-b[2])**2
    )
}
export function closestPointOnBox(ballpos,ballr,boxpos,boxsz){
    const top = boxpos[1]+boxsz[1]/2
    const left = boxpos[0]-boxsz[0]/2
    const right = boxpos[0]+boxsz[0]/2
    const front = boxpos[2]-boxsz[2]/2
    const back = boxpos[2]+boxsz[2]/2
    const bottom = boxpos[1]-boxsz[1]/2
    const corners = {
        ltf: [left,top,front],
        rtf: [right,top,front],
        ltb: [left,top,back],
        rtb: [right,top,back],
        lbf: [left,bottom,front],
        rbf: [right,bottom,front],
        lbb: [left,bottom,back],
        rbb: [right,bottom,back],
    }
    const dists = {
        ltf: vecDist(ballpos, corners.ltf)-ballr,
        rtf: vecDist(ballpos, corners.rtf)-ballr,
        ltb: vecDist(ballpos, corners.ltb)-ballr,
        rtb: vecDist(ballpos, corners.rtb)-ballr,
        lbf: vecDist(ballpos, corners.lbf)-ballr,
        rbf: vecDist(ballpos, corners.rbf)-ballr,
        lbb: vecDist(ballpos, corners.lbb)-ballr,
        rbb: vecDist(ballpos, corners.rbb)-ballr,
    }
    let minDist = Infinity
    let closestPoint = null
    for (const dist in dists){
        if (dist<minDist){
            minDist = dist
            closestPoint = dists[dist]
        }
    }
    // distances from faces
    const faceDists = {
        top: Math.abs(ballpos[1]-top)-ballr,
        bottom: Math.abs(ballpos[1]-bottom)-ballr,
        left: Math.abs(ballpos[0]-left)-ballr,
        right: Math.abs(ballpos[0]-right)-ballr,
        front: Math.abs(ballpos[2]-front)-ballr,
        back: Math.abs(ballpos[2]-back)-ballr,
    }
    for (const dist in faceDists){
        if (faceDists[dist]<minDist){
            minDist = faceDists[dist]
            switch(dist){
                case 'top':
                    closestPoint = [ballpos[0],top,ballpos[2]]
                    break
                case 'bottom':
                    closestPoint = [ballpos[0],bottom,ballpos[2]]
                    break
                case 'left':
                    closestPoint = [left,ballpos[1],ballpos[2]]
                    break
                case 'right':
                    closestPoint = [right,ballpos[1],ballpos[2]]
                    break
                case 'front':
                    closestPoint = [ballpos[0],ballpos[1],front]
                    break
                case 'back':
                    closestPoint = [ballpos[0],ballpos[1],back]
                    break
            }
        }
    }
    return {closestPoint,minDist}
}