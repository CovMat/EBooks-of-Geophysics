BEGIN{
}
{
  if($1=="#"){
    print "#", $15
  }
  else{
    print $0;
  }
}
END{
}
