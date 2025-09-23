
outrd,p

pp=p[0]
az=plasmoid(pp)


device,ret=2
window,1

plmp,p[0]

t=p[0]

t.data=p[1].data^2 + p[2].data^2
t.data=sqrt(t.data)/p[0].data
t.id='P/I'

levels=median(az.data)

levels=levels+ findgen(nint(2*levels))/2

write_png,'i.png',tvrd()

window,2
device,ret=2
plmp,t
;!p.noerase=1
;polvec,p,npts=120
;!p.noerase=0

pp=p[0]
az=plasmoid(pp)

plmp,az,/over,levels=levels,col=cgcolor('red'),/nobar


write_png,'pi.png',tvrd()



atmrd,a

aa=size(a.bz)

bz=a.bz[aa[1]/2,*,*]
by=a.by[aa[1]/2,*,*]

print,minmax(by)
print,minmax(bz)


b=sqrt(bz*bz+by*by)

m=p[4]

m.data=b
window,3
device,ret=2
m.id='b'

plmp,m,title='b'

plmp,az,/over,levels=levels,col=cgcolor('red'),/nobar

write_png,'b.png',tvrd()



blink,[2,3]

end
