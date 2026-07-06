function vol = find_opposite_volume(input_volume,face_projected_ind, ...
    rockmesh)
auxmat = [rockmesh.face.inner.upstream rockmesh.face.inner.downstream; ...
          rockmesh.face.bound.upstream rockmesh.face.bound.upstream];         
auxmat = auxmat(face_projected_ind,:);
ref = ~(auxmat == input_volume);
vol = sum(auxmat .* ref,2);
end
