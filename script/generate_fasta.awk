NR==1 {
    for (i=1; i<=NF; i++) {
        ix[$i] = i
    }
}
NR>1 {
    print ">"$ix[c1]"\n"$ix[c2]
}
