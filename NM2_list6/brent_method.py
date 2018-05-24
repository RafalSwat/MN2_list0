def(x1, xu, tol, f, miter):
    phi =  (1+ np.sqrt(5))/2
    rho = 2 - phi
    u = x1 + rho*(xu - x1)
    v = u
    w = u
    x = u
    fu = f(u)
    fv = fu
    fw = fu
    fx = fu
    xm = 0.5*(x1+xu)
    d = 0 
    e = 0
    for itera in range(1,miter):
        if(abs(x-xm)<=tol):
            break #exit
        else:
            para = abs(e)<tol #zmienna logiczna
            if para:
                r = (x-w)*(fx-fv)
                q = (x-v)*(fx-fw)
                p = (x-v)*q - (x-w)*r
                s = 2*(q-r)
                if(s>0):
                    s = abs(s) #sprawdzamy zachowanie paraboli
                    para = abs(para) < abs(0.5*s*e) and p > s*(x1-x) and p<s*(xu-x)
                    if para:
                        e = d
                        d = p/s
                if not para:
                    if x>= xm:
                        e = x1-x
                    else:
                        e = xu-x
                    d = rho * e
                u = x + d
                fu = f(u)
                if fu <= fx:
                    if u >= x:
                        x1 = x
                    else:
                        xu = x
                    v = w
                    fv = fw
                    w = x
                    fw = fx
                    x = u
                    fx = fu
                else:
                    if u < x:
                        x1 = u
                    else:
                        xu = u
                    
                    if fu <= fw or w==x:
                        v = w
                        fv = fw
                        w = u
                        fw = fu
                    elif fu < fv or v==x or v==w:
                        v = u
                        fv = fu
                xm = 0.5*(x1 + xu)
        output u, fu, itera
