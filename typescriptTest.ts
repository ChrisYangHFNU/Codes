for (let group of DataView.pages("#book").groupBy(p => p.genre)){

    DataView.header(3,group.key);
    DataView.table(["Name","Time Read","Rating"],
              group.rows.sort(k => k.rating,'desc')
              .map(k => [k.file.link,k["time-read"],k.rating]))
}