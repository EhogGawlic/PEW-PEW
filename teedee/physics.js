export function vecDist(a,b){
    return Math.sqrt(
        (a[0]-b[0])**2+
        (a[1]-b[1])**2+
        (a[2]-b[2])**2
    )
}
export function subVec(a,b){
    return [
        a[0]-b[0],
        a[1]-b[1],
        a[2]-b[2]
    ]
}
export function dvc(v,n){
    return [
        v[0]/n,
        v[1]/n,
        v[2]/n
    ]
}
export function mvc(v,n){
    return [
        v[0]*n,
        v[1]*n,
        v[2]*n
    ]
}